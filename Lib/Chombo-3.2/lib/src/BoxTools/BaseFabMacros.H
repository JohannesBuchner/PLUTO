#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// note that DIM>3 versions of these macros require that the "BoxIterator.H"
// is included in the file calling the macros

// Can't use ordinary include guards here or multidim won't work.  So instead
// we do something fancier, with LAST_BASEFABMACROS_H_SPACEDIM:
//#ifndef _BASEFABMACROS_H_
//#define _BASEFABMACROS_H_

#ifndef LAST_BASEFABMACROS_H_SPACEDIM
#define LAST_BASEFABMACROS_H_SPACEDIM 0
#endif
#if CH_SPACEDIM != LAST_BASEFABMACROS_H_SPACEDIM
#undef LAST_BASEFABMACROS_H_SPACEDIM
// #endif is at bottom of file.

#include "SPACE.H"
#if CH_SPACEDIM > 3
#include "BoxIterator.H"
#endif

#ifdef DOXYGEN
#define CH_SPACEDIM  1
#endif

//
// CH_MULTIDIM stuff
//
#ifdef ForAllThisCPencil
#   undef ForAllThisCPencil
#   undef ForAllThisPencil
#   undef ForAllXBPencil
#   undef ForAllXBNNnoindx
#   undef ForAllXBNN
#   undef ForAllXCBNN
#   undef ForAllThisBNN
#   undef ForAllThisCBNN
#   undef ForAllThisBNNXC
#   undef ForAllThisCBNNXC
#   undef ForAllThisBNNXCBN
#   undef ForAllThisBNNXCBNYCBN
#   undef ForAllRevXBNYCBNNN
#   undef EndFor
#   undef EndForTX
#   undef EndForPencil
#endif

#if CH_SPACEDIM == 1
#   define ForAllThisCPencil ForAllThisCPencil1
#   define ForAllThisPencil ForAllThisPencil1
#   define ForAllXBPencil ForAllXBPencil1
#   define ForAllXBNN ForAllXBNN1
#   define ForAllXBNNnoindx ForAllXBNNnoindx1
#   define ForAllXCBNN ForAllXCBNN1
#   define ForAllThisBNN ForAllThisBNN1
#   define ForAllThisCBNN ForAllThisCBNN1
#   define ForAllThisBNNXC ForAllThisBNNXC1
#   define ForAllThisCBNNXC ForAllThisCBNNXC1
#   define ForAllThisBNNXCBN ForAllThisBNNXCBN1
#   define ForAllThisBNNXCBNYCBN ForAllThisBNNXCBNYCBN1
#   define ForAllRevXBNYCBNNN ForAllRevXBNYCBNNN1
#   define EndFor EndFor1
#   define EndForTX EndForTX1
#   define EndForPencil EndForPencil1
#   define LAST_BASEFABMACROS_H_SPACEDIM 1
#endif
#if CH_SPACEDIM == 2
#   define ForAllThisCPencil ForAllThisCPencil2
#   define ForAllThisPencil ForAllThisPencil2
#   define ForAllXBPencil ForAllXBPencil2
#   define ForAllXBNN ForAllXBNN2
#   define ForAllXBNNnoindx ForAllXBNNnoindx2
#   define ForAllXCBNN ForAllXCBNN2
#   define ForAllThisBNN ForAllThisBNN2
#   define ForAllThisCBNN ForAllThisCBNN2
#   define ForAllThisBNNXC ForAllThisBNNXC2
#   define ForAllThisCBNNXC ForAllThisCBNNXC2
#   define ForAllThisBNNXCBN ForAllThisBNNXCBN2
#   define ForAllThisBNNXCBNYCBN ForAllThisBNNXCBNYCBN2
#   define ForAllRevXBNYCBNNN ForAllRevXBNYCBNNN2
#   define EndFor EndFor2
#   define EndForTX EndForTX2
#   define EndForPencil EndForPencil2
#   define LAST_BASEFABMACROS_H_SPACEDIM 2
#endif
#if CH_SPACEDIM == 3
#   define ForAllThisCPencil ForAllThisCPencil3
#   define ForAllThisPencil ForAllThisPencil3
#   define ForAllXBPencil ForAllXBPencil3
#   define ForAllXBNN ForAllXBNN3
#   define ForAllXBNNnoindx ForAllXBNNnoindx3
#   define ForAllXCBNN ForAllXCBNN3
#   define ForAllThisBNN ForAllThisBNN3
#   define ForAllThisCBNN ForAllThisCBNN3
#   define ForAllThisBNNXC ForAllThisBNNXC3
#   define ForAllThisCBNNXC ForAllThisCBNNXC3
#   define ForAllThisBNNXCBN ForAllThisBNNXCBN3
#   define ForAllThisBNNXCBNYCBN ForAllThisBNNXCBNYCBN3
#   define ForAllRevXBNYCBNNN ForAllRevXBNYCBNNN3
#   define EndFor EndFor3
#   define EndForTX EndForTX3
#   define EndForPencil EndForPencil3
#   define LAST_BASEFABMACROS_H_SPACEDIM 3
#endif

#if CH_SPACEDIM > 3
#   define ForAllThisCPencil ForAllThisCPencilHiDim
#   define ForAllThisPencil ForAllThisPencilHiDim
#   define ForAllXBPencil ForAllXBPencilHiDim
#   define ForAllXBNN ForAllXBNNHiDim
#   define ForAllXBNNnoindx ForAllXBNNnoindxHiDim
#   define ForAllXCBNN ForAllXCBNNHiDim
#   define ForAllThisBNN ForAllThisBNNHiDim
#   define ForAllThisCBNN ForAllThisCBNNHiDim
#   define ForAllThisBNNXC ForAllThisBNNXCHiDim
#   define ForAllThisCBNNXC ForAllThisCBNNXCHiDim
#   define ForAllThisBNNXCBN ForAllThisBNNXCBNHiDim
#   define ForAllThisBNNXCBNYCBN ForAllThisBNNXCBNYCBNHiDim
#   define ForAllRevXBNYCBNNN ForAllRevXBNYCBNNNHiDim
#   define EndFor EndForHiDim
#   define EndForTX EndForTXHiDim
#   define EndForPencil EndForPencilHiDim
#   define LAST_BASEFABMACROS_H_SPACEDIM HiDim
#endif

#if (CH_SPACEDIM == 1)

/**
@ingroup macros

  The macro ForAllThisPencil(T,b,ns,nc) is intended to facilitate efficient
  looping over the contents of BaseFabs and objects derived from BaseFab.
  Special attention has been paid to make it work efficiently on vector
  supercomputers.

  This macro acts upon the BaseFab *this.  Instead of iterating over the
  entire Box b, it iterates over components starting at component ns and
  ending at component ns+nc-1, in all directions except the first coordinate
  direction.  The user must iterate over the first coordinate direction within
  the ForAllThisPencil loop.  The macro creates two internal reference
  variables; thisR that references the first element in the pencil, and
  thisLen that gives the length of the pencil.  The first argument of the
  macro is a type: the type contained in the BaseFab that is being iterated
  over.

  We can rewrite the code illustrated in `ForAllThisBNN' in this form:

\code
    template <class T>
    void BaseFab<T>::performSetVal(const T val, const Box bx, int ns, int nc)
    {
       CH_assert(domain.contains(bx));
       ForAllThisPencil(T,bx,ns,nc)
       {
          T* dog = &thisR;
          for (int i = 0; i < thisLen; i++)
             dog[i] = val;
       } EndForPencil
    }
\endcode

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllThisPencil1(T,b,ns,nc)                                    \
  {                                                                     \
    CH_assert(contains(b));                                             \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = this->m_dptr;                                            \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        T *_th_pp = _th_p                                               \
            + ((_b_lo[0] - _th_plo[0])                                  \
               + _n * _th_plen[0]);                                     \
        T &thisR = * _th_pp;                                            \
        const int thisLen = _b_len[0];                                  \

/** @ingroup macros
 Same as ForAllThisPencil, except that it operates on the
 the argument BaseFab passed in. added by bvs,  05/27/99
*/
#define ForAllXBPencil1(T,x,b,ns,nc)                                    \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= x.nComp());                   \
    const int *_th_plo = x.loVect();                                    \
    const IntVect _th_plen = x.size();                                  \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = x.dataPtr();                                             \
    for (int nR = (ns); nR < (ns)+(nc); ++nR)                           \
      {                                                                 \
        T * xR = _th_p                                                  \
            + ((_b_lo[0] - _th_plo[0])                                  \
               + nR * _th_plen[0]);                                     \
    const int thisLen = _b_len[0];

/**  @ingroup macros
  The macro ForAllThisCPencil(T,b,ns,nc) is intended to facilitate efficient
  looping over the contents of BaseFabs and objects derived from BaseFab.
  Special attention has been paid to make it work efficiently on vector
  supercomputers.

  This is the constant version of ForAllThisPencil; i.e. it works when the
  underlying BaseFab is constant.

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllThisCPencil1(T,b,ns,nc)                                   \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _th_p = this->m_dptr;                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        const T *_th_pp = _th_p                                         \
            + ((_b_lo[0] - _th_plo[0])                                  \
               + _n * _th_plen[0]);                                     \
        const T &thisR = * _th_pp;                                      \
        const int thisLen = _b_len[0];

/**  @ingroup macros
  The macro ForAllXBNN(T,x,b,ns,nc) is intended to facilitate efficient
  looping over the contents of BaseFabs and objects derived from BaseFab.
  Special attention has been paid to make it work efficiently on vector
  supercomputers.

  This macro acts upon the BaseFab x where the loop runs over the points in
  the Box b and over components starting at ns and ending at ns+nc-1.  The
  first argument of the macro is a type: the type contained in the BaseFab
  that is being iterated over. The reference variable is xR, where x
  is literally replaced by the macros second argument.  Thus an expression
  ForAllXBNN(int,dog,...) would have a reference variable dogR of type int.

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllXBNN1(T,x,b,ns,nc)                                        \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _x_p = (x) .dataPtr();                                           \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        T *_x_pp = _x_p                                                 \
            + ((_b_lo[0] - _x_plo[0])                                   \
               + _n * _x_plen[0]);                                      \
        for (int _i = 0; _i < _b_len[0]; ++_i, ++_x_pp)                 \
          {                                                             \
            const int iR = _i + _b_lo[0];  (int&)iR += 0;               \
            T &x##R = * _x_pp;

// undocumented interface
// Same as ForALLXBNN() except the index variables are not declared.
// Of course, this assumes they aren't used in the loop.
// This exists only for LevelFluxRegister.cpp, to eliminate the annoying
// 'unused variable' compiler warning msgs.
// dbs  Apr2004
#define ForAllXBNNnoindx1(T,x,b,ns,nc)                                  \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _x_p = (x) .dataPtr();                                           \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        T *_x_pp = _x_p                                                 \
            + ((_b_lo[0] - _x_plo[0])                                   \
               + _n * _x_plen[0]);                                      \
        for (int _i = 0; _i < _b_len[0]; ++_i, ++_x_pp)                 \
          {                                                             \
            T &x##R = * _x_pp;

/**  @ingroup macros
  The macro ForAllXCBNN(T,x,b,ns,nc) is intended to facilitate efficient
  looping over the contents of BaseFabs and objects derived from BaseFab.
  Special attention has been paid to make it work efficiently on vector
  supercomputers.

  This is the constant version of ForAllXBNN; i.e. it works when the
  underlying BaseFab is constant.

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllXCBNN1(T,x,b,ns,nc)                                       \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _x_p = (x).dataPtr();                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n;  (int&)nR += 0;                              \
        const T *_x_pp = _x_p                                           \
            + ((_b_lo[0] - _x_plo[0])                                   \
               + _n * _x_plen[0]);                                      \
        for (int _i = 0; _i < _b_len[0]; ++_i)                          \
          {                                                             \
            const int iR = _i + _b_lo[0]; (int&) iR += 0;               \
            const T & x##R = _x_pp[_i];

/**  @ingroup macros
  The ForAllThisBNN(T,b,ns,nc) macro is intended to facilitate efficient
  looping over the contents of BaseFabs and objects derived from BaseFab.
  Special attention has been paid to make it work efficiently on vector
  supercomputers.

  This macro performs the loop over the current object (*this) where the loop
  runs over the points in the Box b and over components starting at ns and
  ending at ns+nc-1.  The first argument of the macro is a type: the type
  contained in the BaseFab that is being iterated over.  The reference
  variable is thisR.

  For example:
\code
    template<class T>
    void
    BaseFab<T>::performSetVal (const T val, const Box bx, int ns, int num)
    {
      CH_assert(domain.contains(bx));
      ForAllThisBNN(T,bx,ns,num)
      {
        thisR = val;
      } EndFor
    }
\endcode

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllThisBNN1(T,b,ns,nc)                                       \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = this->m_dptr;                                            \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        T *_th_pp = _th_p                                               \
            + ((_b_lo[0] - _th_plo[0])                                  \
               + _n * _th_plen[0]);                                     \
        for (int _i = 0; _i < _b_len[0]; ++_i, ++_th_pp)                \
          {                                                             \
            int iR = _i + _b_lo[0]; iR += 0;                            \
            T &thisR = * _th_pp;

/**  @ingroup macros
  The macro ForAllThisCBNN(T,b,ns,nc) is intended to facilitate efficient
  looping over the contents of BaseFabs and objects derived from BaseFab.
  Special attention has been paid to make it work efficiently on vector
  supercomputers.

  This is the constant version of ForAllThisBNN; i.e. it works when the
  underlying BaseFab is constant.

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllThisCBNN1(T,b,ns,nc)                                      \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _th_p = this->m_dptr;                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        const T *_th_pp = _th_p                                         \
            + ((_b_lo[0] - _th_plo[0])                                  \
               + _n * _th_plen[0]);                                     \
        for (int _i = 0; _i < _b_len[0]; ++_i)                          \
          {                                                             \
            const int iR = _i + _b_lo[0]; (int&)iR += 0;                \
            const T &thisR = _th_pp[_i];

/**  @ingroup macros
  The macro ForAllThisBNNXC(T,b,ns,nc,x,nss) is intended to facilitate
  efficient looping over the contents of BaseFabs and objects derived from
  BaseFab.  Special attention has been paid to make it work efficiently on
  vector supercomputers.

  This macro acts upon the BaseFab *this and in addition is able to utiliize
  values in the const BaseFab x.  The loop runs over the points in the Box b
  and over components starting at ns and ending at ns+nc-1.  The reference
  variables are thisR and xR, respectively.  As usual the x in xR is replaced
  by the macro's fifth argument.  The sixth argument nss is the number of the
  argument in x that corresponds to the ns argument in *this.

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllThisBNNXC1(T,b,ns,nc,x,nss)                               \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            T *_th_pp = _th_p                                           \
                + ((_subbox_lo[0] - _th_plo[0])                         \
                   + _n * _th_plen[0]);                                 \
            const T *_x_pp = _x_p                                       \
                + ((_subbox_lo[0] - _x_plo[0])                          \
                   + _n * _x_plen[0]);                                  \
            for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)       \
              {                                                         \
                int iR = _i + _subbox_lo[0]; iR += 0;                   \
                T &thisR = * _th_pp;                                    \
                const T & x##R = _x_pp[_i];

#define ForAllThisCBNNXC1(T,b,ns,nc,x,nss)                              \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        const T* _th_p = this->dataPtr(ns);                             \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            const T *_th_pp = _th_p                                     \
                + ((_subbox_lo[0] - _th_plo[0])                         \
                   + _n * _th_plen[0]);                                 \
            const T *_x_pp = _x_p                                       \
                + ((_subbox_lo[0] - _x_plo[0])                          \
                   + _n * _x_plen[0]);                                  \
            for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)       \
              {                                                         \
                int iR = _i + _subbox_lo[0]; iR += 0;                   \
                const T &thisR = * _th_pp;                              \
                const T & x##R = _x_pp[_i];

/**  @ingroup macros
  The macro ForAllThisBNNXCBN(T,b,ns,nc,x,bx,nss) is intended to facilitate
  efficient looping over the contents of BaseFabs and objects derived from
  BaseFab.  Special attention has been paid to make it work efficiently on
  vector supercomputers.

  This macro acts upon the BaseFab *this and in addition is able to utiliize
  values in the const BaseFab x.  The loop runs over the points in the
  Box b with components starting at ns and ending at ns+nc-1.  The reference
  variables are thisR and xR, respectively.  As usual the x in xR is replaced
  by the macro's fifth argument.  The sixth argument nss is the number of the
  argument in x that corresponds to the ns argument in *this.  Box bx must
  be the same size as this->box() intersected with b.

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllThisBNNXCBN1(T,b,ns,nc,x,bx,nss)                          \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    CH_assert(bx.sameSize((b)));                                        \
    if (!((b).isEmpty()))                                               \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = (b).loVect();                           \
        const IntVect _subbox_len = (b).size();                         \
        const int *_bx_lo = (bx).loVect();                              \
        T* _th_p = this->dataPtr(ns);                                   \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nss; n##x##R += 0;                       \
            T *_th_pp = _th_p                                           \
                + ((_subbox_lo[0] - _th_plo[0])                         \
                   + _n * _th_plen[0]);                                 \
            const T *_x_pp = _x_p                                       \
                + ((_bx_lo[0] - _x_plo[0])                              \
                   + _n * _x_plen[0]);                                  \
            for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)       \
              {                                                         \
                int iR = _i + _subbox_lo[0]; iR += 0;                   \
                int i##x##R = _i + _bx_lo[0]; i##x##R += 0;             \
                T &thisR = * _th_pp;                                    \
                const T & x##R = _x_pp[_i];

/**  @ingroup macros
  The macro ForAllThisBNNXCBNYCBN(T,b,ns,nc,x,bx,nsx,y,by,nsy) is intended to
  facilitate efficient looping over the contents of BaseFabs and objects
  derived from BaseFab.  Special attention has been paid to make it work
  efficiently on vector supercomputers.

  This macro acts upon the BaseFab *this and in addition is able to utiliize
  values in the const BaseFab x and const BaseFab y.  The loop runs over the
  points in the intersection of Box b with components starting at ns and
  ending at ns+nc-1.  The reference variables are thisR, xR, and yR
  respectively. As usual the x in xR is replaced by the macro's fifth argument
  and likewise for the y in yR.  The seventh argument nsx is the number of the
  argument in x that corresponds to the ns argument in *this, and the eighth
  argument nsy is the number of the argument in y that corresponds to the ns
  argument in *this.  Boxes bx and by must be the same size as this->box()
  intersected with b.

  Looping macro mnemonics:

    This stands for the current object

    C for a const

    X stands for a BaseFab

    B for a Box

    N for an int
*/
#define ForAllThisBNNXCBNYCBN1(T,b,ns,nc,x,bx,nsx,y,by,nsy)             \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    CH_assert((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    Box _subbox_ = this->box();                                         \
    _subbox_ &= b;                                                      \
    CH_assert((bx).sameSize(_subbox_));                                 \
    CH_assert((by).sameSize(_subbox_));                                 \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_y_plo = (y).loVect();                               \
        const IntVect   _y_plen = (y).size();                           \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        const int *_bx_lo = (bx).loVect();                              \
        const int *_by_lo = (by).loVect();                              \
        T* _th_p = this->dataPtr(ns);                                   \
        const T* _x_p  = (x).dataPtr(nsx);                              \
        const T* _y_p  = (y).dataPtr(nsy);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR      = _n + ns;  nR      += 0;                       \
            int n##x##R = _n + nsx; n##x##R += 0;                       \
            int n##y##R = _n + nsy; n##y##R += 0;                       \
            T *_th_pp = _th_p                                           \
                + ((_subbox_lo[0] - _th_plo[0])                         \
                   + _n * _th_plen[0]);                                 \
            const T *_x_pp = _x_p                                       \
                + ((_bx_lo[0] - _x_plo[0])                              \
                   + _n * _x_plen[0]);                                  \
            const T *_y_pp = _y_p                                       \
                + ((_by_lo[0] - _y_plo[0])                              \
                   + _n * _y_plen[0]);                                  \
            for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)       \
              {                                                         \
                int iR = _i + _subbox_lo[0];  iR += 0;                  \
                int i##x##R = _i + _bx_lo[0]; i##x##R += 0;             \
                int i##y##R = _i + _by_lo[0]; i##y##R += 0;             \
                T &thisR = * _th_pp;                                    \
                const T & x##R = _x_pp[_i];                             \
                const T & y##R = _y_pp[_i];

#define ForAllRevXBNYCBNNN1(T,x,bx,nsx,y,by,nsy,nc,ri)                  \
  {                                                                     \
    CH_assert((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    CH_assert((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    CH_assert((x).contains(bx));                                        \
    CH_assert((y).contains(by));                                        \
    CH_assert((bx).sameSize(by));                                       \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_y_plo = (y).loVect();                                   \
    const IntVect _y_plen = (y).size();                                 \
    const IntVect _len = (bx).size();                                   \
    const int *_bx_lo = (bx).loVect();                                  \
    const int *_by_lo = (by).loVect();                                  \
    T* _x_p  = (x).dataPtr(nsx);                                        \
    const T* _y_p  = (y).dataPtr(nsy);                                  \
    for (int _n = 0; _n < (nc); ++_n)                                   \
      {                                                                 \
        int n##x##R = _n + nsx; n##x##R += 0;                           \
        int n##y##R = _n + nsy; n##y##R += 0;                           \
        int _ix = 0;                                                    \
        T *_x_pp = _x_p                                                 \
            + ((_bx_lo[0] - _x_plo[0]) + _len[0] - 1                    \
                + _n * _x_plen[0]);                                     \
        const T *_y_pp = _y_p                                           \
            + ((_by_lo[0] - _y_plo[0])                                  \
                + _n * _y_plen[0]);                                     \
        for (int _i = 0; _i < _len[0]; ++_i, --_ix)                     \
          {                                                             \
            T & x##R = _x_pp[_ix];                                      \
            const T & y##R = _y_pp[_i];

/**  @ingroup macros
  The macro EndForTX must be used to end all ForAllThisBNNXC,
  ForAllThisBNNXCBN and ForAllThisBNNXCBNYCBN looping constructs.
*/
#define EndForTX1 \
              }   \
          }       \
      }           \
  }

/**  @ingroup macros
  The macro EndFor must be used to end all ForAllXBNN, ForAllXCBNN,
  ForAllThisBNN, and ForAllThisCBNN looping constructs.
*/
#define EndFor1 \
          }     \
      }         \
  }

/**  @ingroup macros
  The macro EndForPencil must be used to end ForAll*Pencil looping constructs.
*/
#define EndForPencil1 \
      }               \
  }

#elif (CH_SPACEDIM == 2)

#define ForAllThisCPencil2(T,b,ns,nc)                                   \
  {                                                                     \
    CH_assert(contains(b));                                             \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _th_p = this->m_dptr;                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1];                               \
            const T *_th_pp = _th_p                                     \
                + ((_b_lo[0] - _th_plo[0])                              \
                   + _th_plen[0]*(                                      \
                       (jR - _th_plo[1])                                \
                       + _n * _th_plen[1]));                            \
            const T &thisR = * _th_pp;                                  \
            const int thisLen = _b_len[0];

#define ForAllThisPencil2(T,b,ns,nc)                                    \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = this->m_dptr;                                            \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1]; (int&)jR += 0;                \
            T *_th_pp = _th_p                                           \
                + ((_b_lo[0] - _th_plo[0])                              \
                   + _th_plen[0]*(                                      \
                       (jR - _th_plo[1])                                \
                       + _n * _th_plen[1]));                            \
            T &thisR = * _th_pp;                                        \
            const int thisLen = _b_len[0];                              \

#define ForAllXBPencil2(T,x,b,ns,nc)                                    \
  {                                                                     \
    CH_assert((x).contains(b));                                         \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_th_plo = (x).loVect();                                  \
    const IntVect _th_plen =(x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = (x) .dataPtr();                                          \
    for (int nR = (ns); nR < (ns)+(nc); ++nR)                           \
      {                                                                 \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1];  (int&)jR += 0;               \
            T *xR = _th_p                                               \
                + ((_b_lo[0] - _th_plo[0])                              \
                   + _th_plen[0]*(                                      \
                       (jR - _th_plo[1])                                \
                       + nR * _th_plen[1]));                            \
            const int thisLen = _b_len[0];                              \

#define ForAllXBNN2(T,x,b,ns,nc)                                        \
  {                                                                     \
    CH_assert((x).contains(b));                                         \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _x_p = (x) .dataPtr();                                           \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1]; (int&)jR += 0;                \
            T *_x_pp = _x_p                                             \
                + ((_b_lo[0] - _x_plo[0])                               \
                       + _x_plen[0]*(                                   \
                           (jR - _x_plo[1])                             \
                           + _n * _x_plen[1]));                         \
            for (int _i = 0; _i < _b_len[0]; ++_i, ++_x_pp)             \
              {                                                         \
                const int iR = _i + _b_lo[0];  (int&)iR += 0;           \
                T &x##R = * _x_pp;

#define ForAllXBNNnoindx2(T,x,b,ns,nc)                                  \
  {                                                                     \
    CH_assert((x).contains(b));                                         \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _x_p = (x) .dataPtr();                                           \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1];                               \
            T *_x_pp = _x_p                                             \
                + ((_b_lo[0] - _x_plo[0])                               \
                       + _x_plen[0]*(                                   \
                           (jR - _x_plo[1])                             \
                           + _n * _x_plen[1]));                         \
            for (int _i = 0; _i < _b_len[0]; ++_i, ++_x_pp)             \
              {                                                         \
                T &x##R = * _x_pp;

#define ForAllXCBNN2(T,x,b,ns,nc)                                       \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _x_p = (x).dataPtr();                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1];                               \
            const T *_x_pp = _x_p                                       \
                + ((_b_lo[0] - _x_plo[0])                               \
                       + _x_plen[0]*(                                   \
                           (jR  - _x_plo[1])                            \
                           + _n * _x_plen[1]));                         \
            for (int _i = 0; _i < _b_len[0]; ++_i)                      \
              {                                                         \
                const int iR = _i + _b_lo[0];   (int&)iR += 0;          \
                const T & x##R = _x_pp[_i];

#define ForAllThisBNN2(T,b,ns,nc)                                       \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = this->m_dptr;                                            \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1];                               \
            T *_th_pp = _th_p                                           \
                + ((_b_lo[0] - _th_plo[0])                              \
                   + _th_plen[0]*(                                      \
                       (jR - _th_plo[1])                                \
                       + _n * _th_plen[1]));                            \
            for (int _i = 0; _i < _b_len[0]; ++_i, ++_th_pp)            \
              {                                                         \
                int iR = _i + _b_lo[0]; iR += 0;                        \
                T &thisR = * _th_pp;

#define ForAllThisCBNN2(T,b,ns,nc)                                      \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _th_p = this->m_dptr;                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _j = 0; _j < _b_len[1]; ++_j)                          \
          {                                                             \
            const int jR = _j + _b_lo[1]; (int&)jR+=0;;                 \
            const T *_th_pp = _th_p                                     \
                + ((_b_lo[0] - _th_plo[0])                              \
                   + _th_plen[0]*(                                      \
                       (_j + _b_lo[1] - _th_plo[1])                     \
                       + _n * _th_plen[1]));                            \
            for (int _i = 0; _i < _b_len[0]; ++_i)                      \
              {                                                         \
                const int iR = _i + _b_lo[0]; (int&)iR +=0;             \
                const T &thisR = _th_pp[_i];

#define ForAllThisBNNXC2(T,b,ns,nc,x,nss)                               \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        T* _th_p = this->dataPtr(ns);                                   \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            for (int _j = 0; _j < _subbox_len[1]; ++_j)                 \
              {                                                         \
                const int jR = _j + _subbox_lo[1];                      \
                T *_th_pp = _th_p                                       \
                    + ((_subbox_lo[0] - _th_plo[0])                     \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _n * _th_plen[1]));                        \
                const T *_x_pp = _x_p                                   \
                    + ((_subbox_lo[0] - _x_plo[0])                      \
                       + _x_plen[0]*(                                   \
                           (jR - _x_plo[1])                             \
                           + _n * _x_plen[1]));                         \
                for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)   \
                  {                                                     \
                    int iR = _i + _subbox_lo[0]; iR += 0;               \
                    T &thisR = * _th_pp;                                \
                    const T & x##R = _x_pp[_i];

#define ForAllThisCBNNXC2(T,b,ns,nc,x,nss)                              \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        const T* _th_p = this->dataPtr(ns);                             \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            for (int _j = 0; _j < _subbox_len[1]; ++_j)                 \
              {                                                         \
                const int jR = _j + _subbox_lo[1];                      \
                const T *_th_pp = _th_p                                 \
                    + ((_subbox_lo[0] - _th_plo[0])                     \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _n * _th_plen[1]));                        \
                const T *_x_pp = _x_p                                   \
                    + ((_subbox_lo[0] - _x_plo[0])                      \
                       + _x_plen[0]*(                                   \
                           (jR - _x_plo[1])                             \
                           + _n * _x_plen[1]));                         \
                for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)   \
                  {                                                     \
                    int iR = _i + _subbox_lo[0]; iR += 0;               \
                    const T &thisR = * _th_pp;                          \
                    const T & x##R = _x_pp[_i];

#define ForAllThisBNNXCBN2(T,b,ns,nc,x,bx,nss)                          \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    CH_assert(bx.sameSize((b)));                                        \
    if (!((b).isEmpty()))                                               \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = (b).loVect();                           \
        const IntVect _subbox_len = (b).size();                         \
        const int *_bx_lo = (bx).loVect();                              \
        T* _th_p = this->dataPtr(ns);                                   \
        int nR = ns; int n##x##R = nss;                                 \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n, ++nR, ++n##x##R )             \
          {                                                             \
            for (int _j = 0; _j < _subbox_len[1]; ++_j)                 \
              {                                                         \
                const int jR = _j + _subbox_lo[1];                      \
                const int j##x##R = _j + _bx_lo[1];                     \
                T *_th_pp = _th_p                                       \
                    + ((_subbox_lo[0] - _th_plo[0])                     \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _n * _th_plen[1]));                        \
                const T *_x_pp = _x_p                                   \
                    + ((_bx_lo[0] - _x_plo[0])                          \
                       + _x_plen[0]*(                                   \
                           (j##x##R - _x_plo[1])                        \
                           + _n * _x_plen[1]));                         \
                int iR = _subbox_lo[0]; int i##x##R = _bx_lo[0];        \
                for (int _i = 0; _i < _subbox_len[0];                   \
                     ++_i, ++_th_pp, ++_x_pp, ++iR, ++i##x##R)          \
                  {                                                     \
                    T &thisR = * _th_pp;                                \
                    const T & x##R = * _x_pp;

#define ForAllThisBNNXCBNYCBN2(T,b,ns,nc,x,bx,nsx,y,by,nsy)             \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    CH_assert((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    Box _subbox_ = this->box();                                         \
    _subbox_ &= b;                                                      \
    CH_assert((bx).sameSize(_subbox_));                                 \
    CH_assert((by).sameSize(_subbox_));                                 \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_y_plo = (y).loVect();                               \
        const IntVect _y_plen = (y).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        const int *_bx_lo = (bx).loVect();                              \
        const int *_by_lo = (by).loVect();                              \
        T* _th_p = this->dataPtr(ns);                                   \
        const T* _x_p  = (x).dataPtr(nsx);                              \
        const T* _y_p  = (y).dataPtr(nsy);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nsx; n##x##R += 0;                       \
            int n##y##R = _n + nsy; n##y##R += 0;                       \
            for (int _j = 0; _j < _subbox_len[1]; ++_j)                 \
              {                                                         \
                const int jR = _j + _subbox_lo[1];                      \
                const int j##x##R = _j + _bx_lo[1];                     \
                const int j##y##R = _j + _by_lo[1];                     \
                T *_th_pp = _th_p                                       \
                    + ((_subbox_lo[0] - _th_plo[0])                     \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _n * _th_plen[1]));                        \
                const T *_x_pp = _x_p                                   \
                    + ((_bx_lo[0] - _x_plo[0])                          \
                       + _x_plen[0]*(                                   \
                           (j##x##R - _x_plo[1])                        \
                           + _n * _x_plen[1]));                         \
                const T *_y_pp = _y_p                                   \
                    + ((_by_lo[0] - _y_plo[0])                          \
                       + _y_plen[0]*(                                   \
                           (j##y##R - _y_plo[1])                        \
                           + _n * _y_plen[1]));                         \
                for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)   \
                  {                                                     \
                    int iR = _i + _subbox_lo[0];  iR += 0;              \
                    int i##x##R = _i + _bx_lo[0]; i##x##R += 0;         \
                    int i##y##R = _i + _by_lo[0]; i##y##R += 0;         \
                    T &thisR = * _th_pp;                                \
                    const T & x##R = _x_pp[_i];                         \
                    const T & y##R = _y_pp[_i];

#define ForAllRevXBNYCBNNN2(T,x,bx,nsx,y,by,nsy,nc,ir)                  \
  {                                                                     \
    CH_assert((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    CH_assert((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    CH_assert((ir) >= 0 && (ir) < SpaceDim);                            \
    CH_assert((x).contains(bx));                                        \
    CH_assert((y).contains(by));                                        \
    CH_assert((bx).sameSize(by));                                       \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_y_plo = (y).loVect();                                   \
    const int *_y_plen = (y).size();                                    \
    const int *_bx_lo = (bx).loVect();                                  \
    const int *_by_lo = (by).loVect();                                  \
    const IntVect _len = (bx).size();                                   \
    T* _x_p  = (x).dataPtr(nsx);                                        \
    const T* _y_p  = (y).dataPtr(nsy);                                  \
    for (int _n = 0; _n < (nc); ++_n)                                   \
      {                                                                 \
        int n##x##R = _n + nsx; n##x##R += 0;                           \
        int n##y##R = _n + nsy; n##y##R += 0;                           \
        for (int _j = 0; _j < _len[1]; ++_j)                            \
          {                                                             \
            const int j##x##R = _j + _bx_lo[1];                         \
            const int jrev##x##R = _len[1]-1-_j + _bx_lo[1];            \
            const int j##y##R = _j + _by_lo[1];                         \
            T *_x_pp;                                                   \
            int _ix = 0;                                                \
            int _istrd;                                                 \
            if (ir == 0)                                                \
              {                                                         \
                _x_pp = _x_p                                            \
                    + ((_bx_lo[0] - _x_plo[0]) + _len[0] - 1            \
                       + _x_plen[0]*(                                   \
                           (j##x##R - _x_plo[1])                        \
                           + _n * _x_plen[1]));                         \
                _istrd = -1;                                            \
              }                                                         \
            else                                                        \
              {                                                         \
                _x_pp = _x_p                                            \
                    + ((_bx_lo[0] - _x_plo[0])                          \
                       + _x_plen[0]*(                                   \
                           (jrev##x##R - _x_plo[1])                     \
                           + _n * _x_plen[1]));                         \
                _istrd = 1;                                             \
              }                                                         \
            const T *_y_pp = _y_p                                       \
                    + ((_by_lo[0] - _y_plo[0])                          \
                       + _y_plen[0]*(                                   \
                           (j##y##R - _y_plo[1])                        \
                           + _n * _y_plen[1]));                         \
            int _x_rev = _len[0]-1; _x_rev += 0;                        \
            for (int _i = 0; _i < _len[0]; ++_i, _ix+=_istrd)           \
              {                                                         \
                T & x##R = _x_pp[_ix];                                  \
                const T & y##R = _y_pp[_i];

#define EndFor2 \
              } \
          }     \
      }         \
  }

#define EndForTX2   \
                  } \
              }     \
          }         \
      }             \
  }

#define EndForPencil2 \
          }           \
      }               \
  }

#elif (CH_SPACEDIM == 3)

#define ForAllThisCPencil3(T,b,ns,nc)                                   \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _th_p = this->m_dptr;                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                const T *_th_pp = _th_p                                 \
                    + ((_b_lo[0] - _th_plo[0])                          \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _th_plen[1]*(                              \
                               (kR - _th_plo[2])                        \
                               + _n * _th_plen[2])));                   \
                const T &thisR = * _th_pp;                              \
                const int thisLen = _b_len[0];

#define ForAllThisPencil3(T,b,ns,nc)                                    \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = this->m_dptr;                                            \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                T *_th_pp = _th_p                                       \
                    + ((_b_lo[0] - _th_plo[0])                          \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _th_plen[1]*(                              \
                               (kR - _th_plo[2])                        \
                               + _n * _th_plen[2])));                   \
                T &thisR = * _th_pp;                                    \
                const int thisLen = _b_len[0];                          \

#define ForAllXBPencil3(T,x,b,ns,nc)                                    \
  {                                                                     \
    CH_assert((x).contains(b));                                         \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_th_plo = (x).loVect();                                  \
    const IntVect _th_plen = (x).size();                                \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = (x) .dataPtr();                                          \
    for (int nR = (ns); nR < (ns)+(nc); ++nR)                           \
      {                                                                 \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                T *xR = _th_p                                           \
                    + ((_b_lo[0] - _th_plo[0])                          \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _th_plen[1]*(                              \
                               (kR - _th_plo[2])                        \
                               + nR * _th_plen[2])));                   \
                const int thisLen = _b_len[0];

#define ForAllXBNN3(T,x,b,ns,nc)                                        \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _x_p = (x) .dataPtr();                                           \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                T *_x_pp = _x_p                                         \
                    + ((_b_lo[0] - _x_plo[0])                           \
                       + _x_plen[0]*(                                   \
                           (jR - _x_plo[1])                             \
                           + _x_plen[1]*(                               \
                               (kR - _x_plo[2])                         \
                               + _n * _x_plen[2])));                    \
                for (int _i = 0; _i < _b_len[0]; ++_i, ++_x_pp)         \
                  {                                                     \
                    const int iR = _i + _b_lo[0];  (int&)iR += 0;       \
                    T &x##R = * _x_pp;

#define ForAllXBNNnoindx3(T,x,b,ns,nc)                                  \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _x_p = (x) .dataPtr();                                           \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                T *_x_pp = _x_p                                         \
                    + ((_b_lo[0] - _x_plo[0])                           \
                       + _x_plen[0]*(                                   \
                           (jR - _x_plo[1])                             \
                           + _x_plen[1]*(                               \
                               (kR - _x_plo[2])                         \
                               + _n * _x_plen[2])));                    \
                for (int _i = 0; _i < _b_len[0]; ++_i, ++_x_pp)         \
                  {                                                     \
                    T &x##R = * _x_pp;

#define ForAllXCBNN3(T,x,b,ns,nc)                                       \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _x_p = (x).dataPtr();                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                const T *_x_pp = _x_p                                   \
                    + ((_b_lo[0] - _x_plo[0])                           \
                       + _x_plen[0]*(                                   \
                           (jR  - _x_plo[1])                            \
                           + _x_plen[1]*(                               \
                               (kR - _x_plo[2])                         \
                               + _n * _x_plen[2])));                    \
                for (int _i = 0; _i < _b_len[0]; ++_i)                  \
                  {                                                     \
                    const int iR = _i + _b_lo[0]; (int&)iR += 0;        \
                    const T & x##R = _x_pp[_i];

#define ForAllThisBNN3(T,b,ns,nc)                                       \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    T* _th_p = this->m_dptr;                                            \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                T *_th_pp = _th_p                                       \
                    + ((_b_lo[0] - _th_plo[0])                          \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _th_plen[1]*(                              \
                               (kR - _th_plo[2])                        \
                               + _n * _th_plen[2])));                   \
                for (int _i = 0; _i < _b_len[0]; ++_i, ++_th_pp)        \
                  {                                                     \
                    int iR = _i + _b_lo[0]; iR += 0;                    \
                    T &thisR = * _th_pp;

#define ForAllThisCBNN3(T,b,ns,nc)                                      \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const int *_th_plo = this->loVect();                                \
    const IntVect _th_plen = this->size();                              \
    const int *_b_lo = (b).loVect();                                    \
    const IntVect _b_len = (b).size();                                  \
    const T* _th_p = this->m_dptr;                                      \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        for (int _k = 0; _k < _b_len[2]; ++_k)                          \
          {                                                             \
            const int kR = _k + _b_lo[2];                               \
            for (int _j = 0; _j < _b_len[1]; ++_j)                      \
              {                                                         \
                const int jR = _j + _b_lo[1];                           \
                const T *_th_pp = _th_p                                 \
                    + ((_b_lo[0] - _th_plo[0])                          \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _th_plen[1]*(                              \
                               (kR - _th_plo[2])                        \
                               + _n * _th_plen[2])));                   \
                for (int _i = 0; _i < _b_len[0]; ++_i)                  \
                  {                                                     \
                    const int iR = _i + _b_lo[0]; (int&)iR += 0;        \
                    const T &thisR = _th_pp[_i];

#define ForAllThisBNNXC3(T,b,ns,nc,x,nss)                               \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        T* _th_p = this->dataPtr(ns);                                   \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            for (int _k = 0; _k < _subbox_len[2]; ++_k)                 \
              {                                                         \
                const int kR = _k + _subbox_lo[2];                      \
                for (int _j = 0; _j < _subbox_len[1]; ++_j)             \
                  {                                                     \
                    const int jR = _j + _subbox_lo[1];                  \
                    T *_th_pp = _th_p                                   \
                        + ((_subbox_lo[0] - _th_plo[0])                 \
                           + _th_plen[0]*(                              \
                               (jR - _th_plo[1])                        \
                               + _th_plen[1]*(                          \
                                   (kR - _th_plo[2])                    \
                                   + _n * _th_plen[2])));               \
                    const T *_x_pp = _x_p                               \
                        + ((_subbox_lo[0] - _x_plo[0])                  \
                           + _x_plen[0]*(                               \
                               (jR - _x_plo[1])                         \
                               + _x_plen[1]*(                           \
                                   (kR - _x_plo[2])                     \
                                   + _n * _x_plen[2])));                \
                    for (int _i = 0; _i < _subbox_len[0];               \
                         ++_i, ++_th_pp)                                \
                      {                                                 \
                        int iR = _i + _subbox_lo[0]; iR += 0;           \
                        T &thisR = * _th_pp;                            \
                        const T & x##R = _x_pp[_i];

#define ForAllThisCBNNXC3(T,b,ns,nc,x,nss)                              \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        const T* _th_p = this->dataPtr(ns);                             \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            for (int _k = 0; _k < _subbox_len[2]; ++_k)                 \
              {                                                         \
                const int kR = _k + _subbox_lo[2];                      \
                for (int _j = 0; _j < _subbox_len[1]; ++_j)             \
                  {                                                     \
                    const int jR = _j + _subbox_lo[1];                  \
                    const T *_th_pp = _th_p                             \
                        + ((_subbox_lo[0] - _th_plo[0])                 \
                           + _th_plen[0]*(                              \
                               (jR - _th_plo[1])                        \
                               + _th_plen[1]*(                          \
                                   (kR - _th_plo[2])                    \
                                   + _n * _th_plen[2])));               \
                    const T *_x_pp = _x_p                               \
                        + ((_subbox_lo[0] - _x_plo[0])                  \
                           + _x_plen[0]*(                               \
                               (jR - _x_plo[1])                         \
                               + _x_plen[1]*(                           \
                                   (kR - _x_plo[2])                     \
                                   + _n * _x_plen[2])));                \
                    for (int _i = 0; _i < _subbox_len[0];               \
                         ++_i, ++_th_pp)                                \
                      {                                                 \
                        int iR = _i + _subbox_lo[0]; iR += 0;           \
                        const T &thisR = * _th_pp;                      \
                        const T & x##R = _x_pp[_i];

#define ForAllThisBNNXCBN3(T,b,ns,nc,x,bx,nss)                          \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    CH_assert((bx).sameSize((b)));                                      \
    if (!((b).isEmpty()))                                               \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_subbox_lo = (b).loVect();                           \
        IntVect _subbox_len = (b).size();                               \
        const int *_bx_lo = (bx).loVect();                              \
        T* _th_p = this->dataPtr(ns);                                   \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nss; n##x##R += 0;                       \
            for (int _k = 0; _k < _subbox_len[2]; ++_k)                 \
              {                                                         \
                const int kR = _k + _subbox_lo[2];                      \
                const int k##x##R = _k + _bx_lo[2];                     \
                for (int _j = 0; _j < _subbox_len[1]; ++_j)             \
                  {                                                     \
                    const int jR = _j + _subbox_lo[1];                  \
                    const int j##x##R = _j + _bx_lo[1];                 \
                    T *_th_pp = _th_p                                   \
                        + ((_subbox_lo[0] - _th_plo[0])                 \
                           + _th_plen[0]*(                              \
                               (jR - _th_plo[1])                        \
                               + _th_plen[1]*(                          \
                                   (kR - _th_plo[2])                    \
                                   + _n * _th_plen[2])));               \
                    const T *_x_pp = _x_p                               \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                    for (int _i = 0; _i < _subbox_len[0];               \
                         ++_i, ++_th_pp)                                \
                      {                                                 \
                        int iR = _i + _subbox_lo[0]; iR += 0;           \
                        int i##x##R = _i + _bx_lo[0]; i##x##R += 0;     \
                        T &thisR = * _th_pp;                            \
                        const T & x##R = _x_pp[_i];

#define ForAllThisBNNXCBNYCBN3(T,b,ns,nc,x,bx,nsx,y,by,nsy)             \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    CH_assert((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    Box _subbox_(this->box());                                          \
    _subbox_ &= b;                                                      \
    CH_assert((bx).sameSize(_subbox_));                                 \
    CH_assert((by).sameSize(_subbox_));                                 \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        const int *_th_plo = this->loVect();                            \
        const IntVect _th_plen = this->size();                          \
        const int *_x_plo = (x).loVect();                               \
        const IntVect _x_plen = (x).size();                             \
        const int *_y_plo = (y).loVect();                               \
        const IntVect _y_plen = (y).size();                             \
        const int *_subbox_lo = _subbox_.loVect();                      \
        const IntVect _subbox_len = _subbox_.size();                    \
        const int *_bx_lo = (bx).loVect();                              \
        const int *_by_lo = (by).loVect();                              \
        T* _th_p = this->dataPtr(ns);                                   \
        const T* _x_p  = (x).dataPtr(nsx);                              \
        const T* _y_p  = (y).dataPtr(nsy);                              \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nsx; n##x##R += 0;                       \
            int n##y##R = _n + nsy; n##y##R += 0;                       \
            for (int _k = 0; _k < _subbox_len[2]; ++_k)                 \
              {                                                         \
                const int kR = _k + _subbox_lo[2];                      \
                const int k##x##R = _k + _bx_lo[2];                     \
                const int k##y##R = _k + _by_lo[2];                     \
                for (int _j = 0; _j < _subbox_len[1]; ++_j)             \
                  {                                                     \
                    const int jR = _j + _subbox_lo[1];                  \
                    const int j##x##R = _j + _bx_lo[1];                 \
                    const int j##y##R = _j + _by_lo[1];                 \
                    T *_th_pp = _th_p                                   \
                        + ((_subbox_lo[0] - _th_plo[0])                 \
                           + _th_plen[0]*(                              \
                               (jR - _th_plo[1])                        \
                               + _th_plen[1]*(                          \
                                   (kR - _th_plo[2])                    \
                                   + _n * _th_plen[2])));               \
                    const T *_x_pp = _x_p                               \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                    const T *_y_pp = _y_p                               \
                        + ((_by_lo[0] - _y_plo[0])                      \
                           + _y_plen[0]*(                               \
                               (j##y##R - _y_plo[1])                    \
                               + _y_plen[1]*(                           \
                                   (k##y##R - _y_plo[2])                \
                                   + _n * _y_plen[2])));                \
                    for (int _i = 0; _i < _subbox_len[0];               \
                         ++_i, ++_th_pp)                                \
                      {                                                 \
                        int iR = _i + _subbox_lo[0];  iR += 0;          \
                        int i##x##R = _i + _bx_lo[0]; i##x##R += 0;     \
                        int i##y##R = _i + _by_lo[0]; i##y##R += 0;     \
                        T &thisR = * _th_pp;                            \
                        const T & x##R = _x_pp[_i];                     \
                        const T & y##R = _y_pp[_i];

#define ForAllRevXBNYCBNNN3(T,x,bx,nsx,y,by,nsy,nc,ir)                  \
  {                                                                     \
    CH_assert((ir) >= 0 && (ir) < SpaceDim);                            \
    CH_assert((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    CH_assert((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    CH_assert((x).contains(bx));                                        \
    CH_assert((y).contains(by));                                        \
    CH_assert((bx).sameSize(by));                                       \
    const int *_x_plo = (x).loVect();                                   \
    const IntVect _x_plen = (x).size();                                 \
    const int *_y_plo = (y).loVect();                                   \
    intVect _y_plen = (y).size();                                       \
    const int *_bx_lo = (bx).loVect();                                  \
    const int *_by_lo = (by).loVect();                                  \
    const IntVect _len = (bx).size();                                   \
    T* _x_p  = (x).dataPtr(nsx);                                        \
    const T* _y_p  = (y).dataPtr(nsy);                                  \
    for (int _n = 0; _n < (nc); ++_n)                                   \
      {                                                                 \
        int n##x##R = _n + nsx; n##x##R += 0;                           \
        int n##y##R = _n + nsy; n##y##R += 0;                           \
        for (int _k = 0; _k < _len[2]; ++_k)                            \
          {                                                             \
            const int k##x##R = _k + _bx_lo[2];                         \
            const int krev##x##R = _len[2]-1-_k + _bx_lo[2];            \
            const int k##y##R = _k + _by_lo[2];                         \
            for (int _j = 0; _j < _len[1]; ++_j)                        \
              {                                                         \
                const int j##x##R = _j + _bx_lo[1];                     \
                const int jrev##x##R = _len[1]-1-_j + _bx_lo[1];        \
                const int j##y##R = _j + _by_lo[1];                     \
                T *_x_pp;                                               \
                int _ix = 0;                                            \
                int _istrd = 1;                                         \
                if (ir == 0)                                            \
                  {                                                     \
                    _x_pp = _x_p                                        \
                        + ((_bx_lo[0] - _x_plo[0]) + _len[0]-1          \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                    _istrd = -1;                                        \
                  }                                                     \
                else if (ir == 1)                                       \
                  {                                                     \
                    _x_pp = _x_p                                        \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (jrev##x##R - _x_plo[1])                 \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                  }                                                     \
                else                                                    \
                  {                                                     \
                    _x_pp = _x_p                                        \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (krev##x##R - _x_plo[2])             \
                                   + _n * _x_plen[2])));                \
                  }                                                     \
                const T *_y_pp = _y_p                                   \
                    + ((_by_lo[0] - _y_plo[0])                          \
                       + _y_plen[0]*(                                   \
                           (j##y##R - _y_plo[1])                        \
                           + _y_plen[1]*(                               \
                               (k##y##R - _y_plo[2])                    \
                               + _n * _y_plen[2])));                    \
                for (int _i = 0; _i < _len[0]; ++_i, _ix += _istrd)     \
                  {                                                     \
                    T & x##R = _x_pp[_ix];                              \
                    const T & y##R = _y_pp[_i];

#define EndFor3     \
                  } \
              }     \
          }         \
      }             \
  }

#define EndForTX3       \
                      } \
                  }     \
              }         \
          }             \
      }                 \
  }

#define EndForPencil3 \
              }       \
          }           \
      }               \
  }

#elif (CH_SPACEDIM > 3)

#define ForAllThisCPencilHiDim(T,b,ns,nc)                               \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const IntVect _b_len = (b).size();                                  \
    Box _pencil_box(b);                                                 \
    _pencil_box.setBig(0,b.smallEnd(0));                                \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        BoxIterator _bit( _pencil_box );                                \
        for (_bit.begin(); _bit.ok(); ++_bit)                           \
          {                                                             \
            IntVect _iv = _bit();                                       \
            const T &thisR = this->operator()(_iv,_n);                  \
            const int thisLen = _b_len[0];                              \
            int jR = _iv[1]; jR += 0;                                   \
            int kR = _iv[2]; kR += 0;

#define ForAllThisPencilHiDim(T,b,ns,nc)                                \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    const IntVect _b_len = (b).size();                                  \
    Box _pencil_box(b);                                                 \
    _pencil_box.setBig(0,b.smallEnd(0));                                \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n;                                              \
        BoxIterator _bit( _pencil_box );                                \
        for (_bit.begin(); _bit.ok(); ++_bit)                           \
          {                                                             \
            IntVect _iv = _bit();                                       \
            T &thisR = this->operator()(_iv,nR);                        \
            const int thisLen = _b_len[0];                              \
            int jR = _iv[1]; jR += 0;                                   \
            int kR = _iv[2]; kR += 0;

#define ForAllXBPencilHiDim(T,x,b,ns,nc)                                \
  {                                                                     \
    CH_assert((x).contains(b));                                         \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    const IntVect _b_len = (b).size();                                  \
    Box _pencil_box(b);                                                 \
    _pencil_box.setBig(0,b.smallEnd(0));                                \
    for (int nR = (ns); nR < (ns)+(nc); ++nR)                           \
      {                                                                 \
        BoxIterator _bit( _pencil_box );                                \
        for (_bit.begin(); _bit.ok(); ++_bit)                           \
          {                                                             \
            IntVect _iv = _bit();                                       \
            T *xR = &(x).operator()(_iv,nR);                            \
            const int thisLen = _b_len[0];                              \
            int jR = _iv[1]; jR += 0;                                   \
            int kR = _iv[2]; kR += 0;

#define ForAllXBNNHiDim(T,x,b,ns,nc)                                    \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n; (int&)nR += 0;                               \
        BoxIterator _bit( (b) );                                        \
        for (_bit.begin(); _bit.ok(); ++_bit)                           \
          {                                                             \
            IntVect _iv = _bit();                                       \
            T &x##R = (x)(_iv,nR);                                      \
            T *_x_pp = &x##R; _x_pp += 0;                               \
            int iR = _iv[0]; iR += 0;                                   \
            int jR = _iv[1]; jR += 0;                                   \
            int kR = _iv[2]; kR += 0;

#define ForAllXBNNnoindxHiDim(T,x,b,ns,nc)                              \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
       BoxIterator _bit( (b) );                                         \
       for (_bit.begin(); _bit.ok(); ++_bit)                            \
         {                                                              \
          IntVect _iv = _bit();                                         \
          T& x##R = (x)(_iv,_n);                                        \
          int iR = _iv[0]; iR += 0;                                     \
          int jR = _iv[1]; jR += 0;                                     \
          int kR = _iv[2]; kR += 0;

#define ForAllXCBNNHiDim(T,x,b,ns,nc)                                   \
  {                                                                     \
    CH_assert(x.contains(b));                                           \
    CH_assert((ns) >= 0 && (ns) + (nc) <= (x).nComp());                 \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n;                                              \
        BoxIterator _bit( (b) );                                        \
        for (_bit.begin(); _bit.ok(); ++_bit)                           \
          {                                                             \
            IntVect _iv = _bit();                                       \
            const T & x##R = (x)(_iv, nR);                              \
            const T *_x_pp = &x##R; _x_pp += 0;                         \
            const int iR = _iv[0]; (int&)iR += 0;                       \
            const int jR = _iv[1]; (int&)jR += 0;                       \
            const int kR = _iv[2]; (int&)kR += 0;

#define ForAllThisBNNHiDim(T,b,ns,nc)                                   \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        int nR = _n; nR += 0;                                           \
        BoxIterator _bit( (b));                                         \
        for (_bit.begin(); _bit.ok(); ++_bit)                           \
          {                                                             \
           IntVect _iv = _bit();                                        \
           T &thisR = this->operator()(_iv,nR);                         \
           int iR = _iv[0]; iR += 0;                                    \
           int jR = _iv[1]; jR += 0;                                    \
           int kR = _iv[2]; kR += 0;

#define ForAllThisCBNNHiDim(T,b,ns,nc)                                  \
  {                                                                     \
    CH_assert(this->contains(b));                                       \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    for (int _n = (ns); _n < (ns)+(nc); ++_n)                           \
      {                                                                 \
        const int nR = _n;                                              \
        BoxIterator _bit( (b) );                                        \
        for (_bit.begin(); _bit.ok(); ++_bit)                           \
          {                                                             \
           IntVect _iv = _bit();                                        \
           const T &thisR = this->operator()(_iv,nR);                   \
           const int iR = _iv[0];  (int&)iR +=0;                        \
           const int jR = _iv[1];  (int&)jR +=0;                        \
           const int kR = _iv[2];  (int&)kR +=0;

#define ForAllThisBNNXCHiDim(T,b,ns,nc,x,nss)                           \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            BoxIterator _bit(_subbox_);                                 \
            for (_bit.begin(); _bit.ok(); ++_bit)                       \
              {                                                         \
               IntVect _iv = _bit();                                    \
               T &thisR = this->operator()(_iv, ns+nR);                 \
               const T & x##R = (x)(_iv, nss+nR);                       \
               const int iR = _iv[0]; (int&)iR += 0;                    \
               const int jR = _iv[1]; (int&)jR += 0;                    \
               const int kR = _iv[2]; (int&)kR += 0;

#define ForAllThisCBNNXCHiDim(T,b,ns,nc,x,nss)                          \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= this->box();                                            \
    _subbox_ &= b;                                                      \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n; nR += 0;                                       \
            BoxIterator _bit(_subbox_);                                 \
            for (_bit.begin(); _bit.ok(); ++_bit)                       \
              {                                                         \
               IntVect _iv = _bit();                                    \
               const T &thisR = this->operator()(_iv,ns+nR);            \
               const T & x##R = (x)(_iv,nss+nR);                        \
               const int iR = _iv[0]; (int&)iR += 0;                    \
               const int jR = _iv[1]; (int&)jR += 0;                    \
               const int kR = _iv[2]; (int&)kR += 0;

#define ForAllThisBNNXCBNHiDim(T,b,ns,nc,x,bx,nss)                      \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    CH_assert((bx).sameSize((b)));                                      \
    if (!((b).isEmpty()))                                               \
      {                                                                 \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nss; n##x##R += 0;                       \
            BoxIterator this_Bit( (b) );                                \
            IntVect _bxOffset = bx.smallEnd() - b.smallEnd();           \
            for (this_Bit.begin(); this_Bit.ok(); ++this_Bit)           \
              {                                                         \
               IntVect _thisIV = this_Bit();                            \
               IntVect _xIV = _thisIV + _bxOffset;                      \
               T &thisR = this->operator()(_thisIV, nR);                \
               const T & x##R = (x)(_xIV,n##x##R);                      \
               int iR = _xIV[0]; iR += 0;                               \
               int jR = _xIV[1]; jR += 0;                               \
               int kR = _xIV[2]; kR += 0;

#define ForAllThisBNNXCBNYCBNHiDim(T,b,ns,nc,x,bx,nsx,y,by,nsy)         \
  {                                                                     \
    CH_assert((ns) >= 0 && (ns) + (nc) <= this->nComp());               \
    CH_assert((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    CH_assert((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    Box _subbox_(this->box());                                          \
    _subbox_ &= b;                                                      \
    CH_assert((bx).sameSize(_subbox_));                                 \
    CH_assert((by).sameSize(_subbox_));                                 \
    if (!_subbox_.isEmpty())                                            \
      {                                                                 \
        for (int _n = 0; _n < (nc); ++_n)                               \
          {                                                             \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nsx; n##x##R += 0;                       \
            int n##y##R = _n + nsy; n##y##R += 0;                       \
            BoxIterator _thisBit(_subbox_);                             \
            IntVect _xOffset = bx.smallEnd() - _subbox_.smallEnd();     \
            IntVect _yOffset = by.smallEnd() - _subbox_.smallEnd();     \
            for (_thisBit.begin(); _thisBit.ok(); ++_thisBit)           \
              {                                                         \
                IntVect _thisIV = _thisBit();                           \
                IntVect _xIV = _thisIV + _xOffset;                      \
                IntVect _yIV = _thisIV + _yOffset;                      \
                T &thisR = this->operator()(_thisIV,nR);                \
                const T & x##R = (x)(_xIV,n##x##R);                     \
                const T & y##R = (y)(_yIV,n##y##R);                     \
                int iR = _thisIV[0]; iR += 0;                           \
                int jR = _thisIV[1]; jR += 0;                           \
                int kR = _thisIV[2]; kR += 0;

#define EndForHiDim \
          }         \
      }             \
  }

#define EndForTXHiDim \
              }       \
          }           \
      }               \
  }

#define EndForPencilHiDim \
          }         \
      }             \
  }

#endif
/**  @ingroup macros
  The Macro ForAllXPencil is a shortened version of ForAllXBPencil
  where the box is defaulted to the entire domain the the
  components run over all components of x
  */
#define ForAllXPencil(T,x)  ForAllXBPencil(T,x,((x).box()),0,((x).nComp()))

/**  @ingroup macros
  The macro ForAllX(T,x) is a shortened form of `ForAllXBNN' where the Box
  defaults to the domain of x and the components run over all the components
  of x.
*/
#define ForAllX(T,x)            ForAllXBNN(T,x,((x).box()),0,((x).nComp()))

/**  @ingroup macros
  The macro ForAllXC(T,x) is the constant form of ForAllX(T,x).
*/
#define ForAllXC(T,x)           ForAllXCBNN(T,x,((x).box()),0,((x).nComp()))

/**  @ingroup macros
  The macro ForAllXB(T,x,b) is a shortened form of `ForAllXBNN'
  where the components run over all the components of x.
*/
#define ForAllXB(T,x,b)         ForAllXBNN(T,x,(b),0,(x).nComp())

/**  @ingroup macros
  The macro ForAllXBC(T,x,b) is the constant form of ForAllXB(T,x,b).
*/
#define ForAllXBC(T,x,b)        ForAllXCBNN(T,x,(b),0,(x).nComp())

/**  @ingroup macros
  The macro ForAllThis(T) is a shortened form of `ForAllThisBNN' where the Box
  defaults to the domain of x and the components run over all the components
  of x.
*/
#define ForAllThis(T)           ForAllThisBNN(T,this->m_domain,0,this->nComp())

/**  @ingroup macros
  The macro ForAllThisC(T) is the constant form of ForAllThis(T).
*/
#define ForAllThisC(T)          ForAllThisCBNN(T,this->m_domain,0,this->nComp())

/**  @ingroup macros
  The macro ForAllThisB(T,b) is a shortened form of `ForAllThisBNN'
  where the components run over all the components of x.
*/
#define ForAllThisB(T,b)        ForAllThisBNN(T,(b),0,this->nComp())

/**  @ingroup macros
  The macro ForAllThisCB(T,b) is the constant form of ForAllThisB(T,b).
*/
#define ForAllThisCB(T,b)       ForAllThisCBNN(T,(b),0,this->nComp())

/**  @ingroup macros
  The macro ForAllThisNN(T,ns,nc) is a shortened form of `ForAllThisBNN'
  where the Box defaults to the domain of *this.
*/
#define ForAllThisNN(T,ns,nc)     ForAllThisBNN(T,this->m_domain,ns,nc)

/**  @ingroup macros
  The macro ForAllThisXC(T,x) is a shortened form of `ForAllThisBNNXC'
  where the Box defaults to the domain of *this and the components run over
  all the components of *this.
*/
#define ForAllThisXC(T,x)       ForAllThisBNNXC(T,this->m_domain,0,this->nComp(),x,0)

#ifdef DOXYGEN
#undef CH_SPACEDIM
#endif /* Doxygen */

#endif // matches CH_SPACEDIM != LAST_BASEFABMACROS_H_SPACEDIM
