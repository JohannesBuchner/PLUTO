#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <climits>
#include "SPACE.H"

#include "MayDay.H"
#include "Misc.H"
#include "Box.H"
#include "parstream.H"
#include "SliceSpec.H"
#include "NamespaceHeader.H"

using std::cout;
using std::cerr;
using std::endl;
using std::ws;
using std::ostream;
using std::istream;

// const Box Box::Empty = Box(IntVect::Unit,
//                            IntVect::Zero);

CH_XDIR::IntVect
convertOldToNew(const IntVect& a_ivOld,
                const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivNew;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivNew[idir] =
        a_sign[idir]*a_ivOld[a_permutation[idir]] + a_translation[idir];
    }
  return ivNew;
}

CH_XDIR::IntVect
convertNewToOld(const IntVect& a_ivNew,
                const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivOld;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivOld[a_permutation[idir]] = a_sign[idir] *
        (a_ivNew[idir] - a_translation[idir]);
    }
  return ivOld;
}

///multiblock stuff.
void
Box::
convertOldToNew(const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivNewLo = CH_XDIR::convertOldToNew(smallEnd(), a_permutation, a_sign, a_translation);
  IntVect ivNewHi = CH_XDIR::convertOldToNew(bigEnd()  , a_permutation, a_sign, a_translation);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int iLo = Min(ivNewLo[idir], ivNewHi[idir]);
      int iHi = Max(ivNewLo[idir], ivNewHi[idir]);
      ivNewLo[idir] = iLo;
      ivNewHi[idir] = iHi;
    }
  //  Box bxNewNodes(ivNewLo, ivNewHi, IndexType::TheNodeType());
  //  Box bxNewCells = enclosedCells(bxNewNodes);
  *this = Box(ivNewLo, ivNewHi);
}

///multiblock stuff
void
Box::
convertNewToOld(const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivOldLo = CH_XDIR::convertNewToOld(smallEnd(), a_permutation, a_sign, a_translation);
  IntVect ivOldHi = CH_XDIR::convertNewToOld(bigEnd()  , a_permutation, a_sign, a_translation);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int iLo = Min(ivOldLo[idir], ivOldHi[idir]);
      int iHi = Max(ivOldLo[idir], ivOldHi[idir]);
      ivOldLo[idir] = iLo;
      ivOldHi[idir] = iHi;
    }
  //  Box bxOldNodes(ivOldLo, ivOldHi, IndexType::TheNodeType());
  //  Box bxOldCells = enclosedCells(bxOldNodes);
  *this =  Box(ivOldLo, ivOldHi);
}

IndexType
IndexType::TheCellType ()
{
    static const IndexType Cell(D_DECL6(IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL));
    return Cell;
}

IndexType
IndexType::TheNodeType ()
{
    static const IndexType Node(D_DECL6(IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE));

    return Node;
}

ostream&
operator<< (ostream&         os,
            const IndexType& it)
{
    os << '('
       << D_TERM6( (it.test(0)?'N':'C'),
                   << ',' << (it.test(1)?'N':'C'),
                   << ',' << (it.test(2)?'N':'C'),
                   << ',' << (it.test(3)?'N':'C'),
                   << ',' << (it.test(4)?'N':'C'),
                   << ',' << (it.test(5)?'N':'C')) << ')' ;
    if (os.fail())
        MayDay::Error("operator<<(ostream&,IndexType&) failed");

    return os;
}

//
// Copied from <Utility.H>
//
#define CH_IGNORE_MAX 100000

istream&
operator>> (istream&   is,
            IndexType& it)
{
  char D_DECL6(t0,t1,t2,
               t3,t4,t5);

  D_EXPR6( is.ignore(CH_IGNORE_MAX, '(') >> t0,
           is.ignore(CH_IGNORE_MAX, ',') >> t1,
           is.ignore(CH_IGNORE_MAX, ',') >> t2,
           is.ignore(CH_IGNORE_MAX, ',') >> t3,
           is.ignore(CH_IGNORE_MAX, ',') >> t4,
           is.ignore(CH_IGNORE_MAX, ',') >> t5);
  is.ignore(CH_IGNORE_MAX, ')');
  D_TERM6(
          CH_assert(t0 == 'C' || t0 == 'N'); t0=='N'?it.set(0):it.unset(0); ,
          CH_assert(t1 == 'C' || t1 == 'N'); t1=='N'?it.set(1):it.unset(1); ,
          CH_assert(t2 == 'C' || t2 == 'N'); t2=='N'?it.set(2):it.unset(2); ,
          CH_assert(t3 == 'C' || t3 == 'N'); t3=='N'?it.set(3):it.unset(3); ,
          CH_assert(t4 == 'C' || t4 == 'N'); t4=='N'?it.set(4):it.unset(4); ,
          CH_assert(t5 == 'C' || t5 == 'N'); t5=='N'?it.set(5):it.unset(5));

  if (is.fail())
        MayDay::Error("operator>>(ostream&,IndexType&) failed");

    return is;
}

// const Box&
// Box::TheUnitBox ()
// {
//     static const Box Unit(IntVect::Zero, IntVect::Unit);
//     return Unit;
// }

//
// Administrative functions.
//
Box::Box ()
    : smallend(IntVect::Unit),
      bigend(IntVect::Zero),
      //      len(IntVect::Zero),
      btype()
{
}

Box::Box (const IntVect& small,
          const int*     vec_len)
    : smallend(small),
      bigend(D_DECL6(small[0]+vec_len[0]-1,
                     small[1]+vec_len[1]-1,
                     small[2]+vec_len[2]-1,
                     small[3]+vec_len[3]-1,
                     small[4]+vec_len[4]-1,
                     small[5]+vec_len[5]-1))
{
  CH_assert(D_TERM6(vec_len[0] >= 0, && vec_len[1] >= 0, && vec_len[2] >= 0,
                    && vec_len[3] >= 0, && vec_len[4] >= 0, && vec_len[5] >= 0));

  //  D_EXPR(len[0] = vec_len[0], len[1] = vec_len[1], len[2] = vec_len[2]);
}

Box::Box (const IntVect&   small,
          const IntVect&   big,
          const IndexType& t)
    : smallend(small),
      bigend(big),
      btype(t)
{
  CH_assert (small <= big);

  computeBoxLen();
}

void Box::define(const IntVect&   small,
                 const IntVect&   big,
                 const IndexType& t)
{
  CH_assert (small <= big);
  smallend = small;
  bigend   = big;
  btype    = t;

  computeBoxLen();
}

Box::Box (const IntVect& small,
          const IntVect& big)
    : smallend(small),
      bigend(big)
{
  CH_assert (small <= big);

  computeBoxLen();
}

void Box::define (const IntVect& small,
                  const IntVect& big)
{
  smallend = small;
  bigend=big;
  CH_assert (small <= big);
  computeBoxLen();
}

Box::Box (const IntVect& small,
          const IntVect& big,
          const IntVect& typ)
    : smallend(small),
      bigend(big),
      btype(typ)
{
  CH_assert (small <= big);
  CH_assert(typ >= IntVect::Zero && typ <= IntVect::Unit);

  computeBoxLen();
}

void Box::define(const IntVect& small,
                 const IntVect& big,
                 const IntVect& typ)
{
  CH_assert (small <= big);
  CH_assert(typ >= IntVect::Zero && typ <= IntVect::Unit);
  smallend = small;
  bigend   = big;
  btype    = typ;
  computeBoxLen();
}

// note that if b is undefined, then this will be undefined as well
// (but is not an error)
void Box::define(const Box& b)
{
  smallend = b.smallend;
  bigend = b.bigend;
  //  len = b.len;
  btype = b.btype;
}

bool
Box::numPtsOK (long& N) const
{

#if CH_SPACEDIM >= 1
    long M = (bigend[0]-smallend[0]+1);
    N = M;
#endif
#if CH_SPACEDIM >= 2
    M = (bigend[1]-smallend[1]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 3
    M = (bigend[2]-smallend[2]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 4
    M = (bigend[3]-smallend[3]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 5
    M = (bigend[4]-smallend[4]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 6
    M = (bigend[5]-smallend[5]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM < 1 || CH_SPACEDIM > 6
    MayDay::Error("Box::numPtsOK() not implmented");
#endif

    return true;
}

long
Box::numPts () const
{
  if (bigend[0] < smallend[0]) return 0;

  long result;
  if (!numPtsOK(result))
    {
      cerr << "Bad Box::numPts:  box = " << *this << endl;
      MayDay::Error("Arithmetic overflow in Box::numPts()");
    }
  return result;
}

bool
Box::volumeOK (long& N) const
{
  return(numPtsOK(N));
}

long
Box::volume () const
{
  return(numPts());
}

Box&
Box::shiftHalf (int dir,
                int nzones)
{
  if (!isEmpty())
  {
    int nbit = (nzones<0 ? -nzones : nzones)%2;
    int nshift = nzones/2;
    unsigned int bit_dir = btype[dir];
    //
    // Toggle btyp bit if nzones is odd.
    //
    if (nbit)
        btype.flip(dir);
    if (nzones < 0)
        nshift -= (bit_dir ? nbit : 0);
    else
        nshift += (bit_dir ? 0 : nbit);
    smallend.shift(dir,nshift);
    bigend.shift(dir,nshift);
  }
  return *this;
}

//
// Boolean functions.
//

bool
Box::intersects (const Box& b) const
{
  if (b.isEmpty() || isEmpty()) return false;
  CH_assert(sameType(b));
  IntVect low(smallend);
  IntVect hi(bigend);
  low.max(b.smallend);
  hi.min(b.bigend);
  return low <= hi;
}

// do boxes intersect?
bool Box::intersectsNotEmpty(const Box& b) const
{
  IntVect low(smallend);
  IntVect hi(bigend);
  low.max(b.smallend);
  hi.min(b.bigend);
  return low <= hi;
}

//
// Intersection functions.
//

Box&
Box::operator&= (const Box& b)
{

  CH_assert(sameType(b));
  if (isEmpty() || b.isEmpty())
  {
    *this = Box().convert(b.btype);
  }
  else
  {
    smallend.max(b.smallend);
    bigend.min(b.bigend);
    if (!(bigend >= smallend)) *this = Box().convert(b.btype);
    else computeBoxLen(); // I'd like to stop doing this....bvs
  }

  return *this;
}

Box&
Box::surroundingNodes ()
{
  if (!isEmpty())
  {
    for (int i = 0; i < CH_SPACEDIM; ++i)
        if ((btype[i] == 0))
            bigend.shift(i,1);
  }
    btype.setall();
    computeBoxLen();
    return *this;
}

Box&
Box::enclosedCells ()
{
  if (!isEmpty())
  {
    for (int i = 0 ; i < CH_SPACEDIM; ++i)
        if (btype[i])
            bigend.shift(i,-1);
  }
    btype.clear();
    computeBoxLen();
    return *this;
}

//
// Next:step through the rectangle with unit increment.
//

void
Box::next (IntVect& p) const
{
    CH_assert(contains(p));

    p.shift(0,1);
#if CH_SPACEDIM == 1
    // nothing more to do
#elif CH_SPACEDIM==2
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
    }
#elif CH_SPACEDIM==3
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
        }
    }
#elif CH_SPACEDIM == 4
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,1);
              }
        }
    }
#elif CH_SPACEDIM == 5
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,1);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,1);
                  }
              }
        }
    }
#elif CH_SPACEDIM == 6
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,1);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,1);
                    if (!(p <= bigend))
                      {
                        p.setVal(4,smallend[4]);
                        p.shift(5,1);
                      }
                  }
              }
        }
    }
#else
    bad spacedim!
#endif
}

//
// Scan point over region of object Box with a vector incrment
// point incrments by 0 direction portion of vector.  When end of
// Box is reached, an increment is made with the 1 direction portion
// of the vector, and the 0 direction scan is resumed.
// effectively, we are scanning a Box, whose length vector is the argument
// vector over the object Box.
// when scan is over, the argument point is over edge of object Box
// this is the signal that we can go no further.
//

void
Box::next (IntVect&   p,
           const int* shv) const
{
    CH_assert(contains(p));

#if   CH_SPACEDIM==1
    p.shift(0,shv[0]);
#elif CH_SPACEDIM==2
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
    }
#elif CH_SPACEDIM==3
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
        }
    }
#elif CH_SPACEDIM==4
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,shv[3]);
              }
        }
    }
#elif CH_SPACEDIM==5
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,shv[3]);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,shv[4]);
                  }
              }
        }
    }
#elif CH_SPACEDIM==6
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,shv[3]);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,shv[4]);
                    if (!(p <= bigend))
                      {
                        p.setVal(4,smallend[4]);
                        p.shift(5,shv[5]);
                      }
                  }
              }
        }
    }
#else
    bad_spacedim!
#endif
}

Box&
Box::refine (int refinement_ratio)
{
  CH_assert(refinement_ratio > 0);
  if (!isEmpty())
  {
    IntVect shft(IntVect::Unit);
    shft -= btype.ixType();
    smallend.scale(refinement_ratio);
    bigend += shft; // Bigend does more than just multiply.
    bigend.scale(refinement_ratio);
    bigend -= shft;
    computeBoxLen();
  }
    return *this;
}

Box&
Box::refine (const IntVect& refinement_ratio)
{
  if (!isEmpty())
  {
    IntVect shft(IntVect::Unit);
    shft -= btype.ixType();
    smallend *= refinement_ratio;
    bigend += shft;
    bigend   *= refinement_ratio;
    bigend -= shft;
    computeBoxLen();
  }
    return *this;
}

//
// Define a macro which will compute an object's length vector from
// the smallend and bigend.  Do several versions according to dimension
// requires you to be in a member functions.

int
Box::longside () const
{
  int maxlen = bigend[0]-smallend[0]+1;
  for (int i = 1; i < CH_SPACEDIM; i++)
    {
      int len=bigend[i]-smallend[i]+1;
        if (len > maxlen)
          maxlen = len;
    }
  return maxlen;
}

int
Box::longside (int& dir) const
{
  int maxlen = bigend[0]-smallend[0]+1;
  dir=0;
  for (int i = 1; i < CH_SPACEDIM; i++)
    {
      int len=bigend[i]-smallend[i]+1;
        if (len > maxlen)
          {
            maxlen = len;
            dir=i;
          }
    }
  return maxlen;
}

int
Box::shortside () const
{
  int dir;
  return shortside(dir);
}

int
Box::shortside (int& dir) const
{
  int minlen = bigend[0]-smallend[0]+1;
  dir=0;
  for (int i = 1; i < CH_SPACEDIM; i++)
    {
      int len=bigend[i]-smallend[i]+1;
        if (len < minlen)
          {
            minlen = len;
            dir=i;
          }
    }
  return minlen;
}

//
// Modified Box is low end, returned Box is high end.
// If CELL: chop_pnt included in high end.
// If NODE: chop_pnt included in both Boxes.
//

Box
Box::chop (int dir,
           int chop_pnt)
{
  CH_assert(!isEmpty());
  //
  // Define new high end Box including chop_pnt.
  //
  IntVect sm(smallend);
  IntVect bg(bigend);
  sm.setVal(dir,chop_pnt);
  if (btype[dir])
  {
    //
    // NODE centered Box.
    //
    CH_assert(chop_pnt > smallend[dir] && chop_pnt < bigend[dir]);
    //
    // Shrink original Box to just contain chop_pnt.
    //
    bigend.setVal(dir,chop_pnt);
  }
  else
  {
    //
    // CELL centered Box.
    //
    CH_assert(chop_pnt > smallend[dir] && chop_pnt <= bigend[dir]);
    //
    // Shrink origional Box to one below chop_pnt.
    //
    bigend.setVal(dir,chop_pnt-1);
  }
  computeBoxLen();
  return Box(sm,bg,btype);
}

void
Box::degenerate( Box& a_to, const SliceSpec& a_sliceSpec,
                 bool* a_outofbounds /*=0*/ ) const
{
  CH_assert(!isEmpty());
  CH_assert( (a_sliceSpec.direction >= 0) && (a_sliceSpec.direction < CH_SPACEDIM) );
  IntVect smTo, bgTo;
  for ( int i=0;i<CH_SPACEDIM;++i )
    {
      smTo[i] = this->smallend[i];
      bgTo[i] = this->bigend[i];
    }
  smTo[a_sliceSpec.direction] = a_sliceSpec.position;
  bgTo[a_sliceSpec.direction] = a_sliceSpec.position;
  Box result( smTo, bgTo );

  if ( a_outofbounds )
    {
      if ( (a_sliceSpec.position < this->smallend[a_sliceSpec.direction])
      ||  (a_sliceSpec.position > this->bigend[a_sliceSpec.direction]) )
        {
          *a_outofbounds = true;
        }
      else
        {
          *a_outofbounds = false;
        }
    }
  a_to = result;
}

//
// Shift by half increments.
//

Box&
Box::shiftHalf (const IntVect& nz)
{
  if (!isEmpty())
  {
   for (int dir = 0; dir < CH_SPACEDIM; dir++)
   {
       int nzones = nz[dir];
       int nbit = (nzones<0 ? -nzones : nzones)%2;
       int nshift = nzones/2;
       unsigned int bit_dir = btype[dir];
       //
       // Toggle btype bit if nzones is odd.
       //
       if (nbit)
           btype.flip(dir);
       if (nzones < 0)
           nshift -= (bit_dir ? nbit : 0);
       else
           nshift += (bit_dir ? 0 : nbit);
       smallend.shift(dir,nshift);
       bigend.shift(dir,nshift);
   }
  }
   return *this;
}

Box& Box::shift_intvect (const IntVect& iv)
{
  return shift(iv);
}

Box& Box::enclosedCells_int (int dir)
{
  return enclosedCells(dir);
}

Box& Box::surroundingNodes_int(int dir)
{
  return surroundingNodes(dir);
}

bool Box::lt(const Box& rhs) const
{
  return *this < rhs;
}

bool Box::eq(const Box& b) const
{
  return *this == b;
}
bool Box::neq(const Box& b) const
{
  return *this != b;
}

Box& Box::shiftHalf_intvect (const IntVect& iv)
{
  return shiftHalf(iv);
}

Box&
Box::convert (int                  dir,
              IndexType::CellIndex typ)
{
  if (!isEmpty())
  {
   unsigned int bitval = btype[dir];
   int off = typ - bitval;
   bigend.shift(dir,off);
   if (off != 0)
      computeBoxLen();
  }
   //
   // Set dir'th bit to typ.
   //
   btype.setType(dir,typ);
   return *this;
}

Box&
Box::convert (IndexType t)
{
   for (int dir = 0; dir < CH_SPACEDIM; dir++)
   {
      unsigned int typ = t[dir];
      if (!isEmpty())
      {
        unsigned int bitval = btype[dir];
        int off = typ - bitval;
        bigend.shift(dir,off);
      }
      btype.setType(dir, (IndexType::CellIndex) typ);
   }
   computeBoxLen();
   return *this;
}

//
// Refinement functions.
//

Box
refine (const Box& b,
        int        refinement_ratio)
{
  CH_assert(refinement_ratio > 0);
  if (b.isEmpty()) return(b);

    IntVect small(b.smallend);
    small.scale(refinement_ratio);
    IntVect big(b.bigend);
    IntVect shft(IntVect::Unit);
    shft -= b.btype.ixType();
    big += shft;        // Large end is not just multiply.
    big.scale(refinement_ratio);
    big -= shft;
    return Box(small,big,b.btype);
}

Box
refine (const Box&     b,
        const IntVect& refinement_ratio)
{
  if (b.isEmpty()) return(b);

    IntVect small(b.smallend);
    small *= refinement_ratio;
    IntVect big(b.bigend);
    IntVect shft(IntVect::Unit);
    shft -= b.btype.ixType();
    big += shft;
    big *= refinement_ratio;
    big -= shft;
    return Box(small,big,b.btype);
}

//
// Coarsening.
//

Box&
Box::coarsen (const IntVect& refinement_ratio)
{
  if (!isEmpty())
  {
    smallend.coarsen(refinement_ratio);

    if (btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (btype[dir])
            {
                int b = bigend[dir];
                int r = refinement_ratio[dir];
                if (b%r)
                    off.setVal(dir,1);
            }
        }
        bigend.coarsen(refinement_ratio);
        bigend += off;
    }
    else
        bigend.coarsen(refinement_ratio);

    computeBoxLen();
  }
    return *this;
}

Box
coarsen (const Box& b,
         int  refinement_ratio)
{
  if (b.isEmpty()) return(b);

    IntVect small(b.smallend);
    small.coarsen(refinement_ratio);
    IntVect big(b.bigend);

    if (b.btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (b.btype[dir])
                if (big[dir]%refinement_ratio)
                    off.setVal(dir,1);
        }
        big.coarsen(refinement_ratio);
        big += off;
    }
    else
        big.coarsen(refinement_ratio);

    return Box(small,big,b.btype);
}

Box&
Box::coarsen (int refinement_ratio)
{
  if (!isEmpty())
  {
    smallend.coarsen(refinement_ratio);
    if (btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (btype[dir])
                if (bigend[dir]%refinement_ratio)
                    off.setVal(dir,1);
        }
        bigend.coarsen(refinement_ratio);
        bigend += off;
    }
    else
        bigend.coarsen(refinement_ratio);

    computeBoxLen();
  }
    return *this;
}

Box
coarsen (const Box&     b,
         const IntVect& refinement_ratio)
{
  if (b.isEmpty()) return(b);

    IntVect small(b.smallend);
    small.coarsen(refinement_ratio);
    IntVect big(b.bigend);

    if (b.btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (b.btype[dir])
                if (big[dir]%refinement_ratio[dir])
                    off.setVal(dir,1);
        }
        big.coarsen(refinement_ratio);
        big += off;
    }
    else
        big.coarsen(refinement_ratio);

    return Box(small,big,b.btype);
}

//
// I/O functions.
//

ostream&
operator<< (ostream&   os,
            const Box& b)
{
  if ( Box::s_tempestOutputFormat )
  {
    os << '['
       << b.smallend << ','
       << b.bigend
       << ']';
  } else
  {
    os << '('
       << b.smallend << ' '
       << b.bigend   << ' '
       << b.btype.ixType()
       << ')';
  }
  if (os.fail())
      MayDay::Error("operator<<(ostream&,Box&) failed");
  return os;
}

void Box::p() const
{
  pout() << *this << "\n";
}

//
// Moved out of Utility.H
//
#define CH_IGNORE_MAX 100000

istream&
operator>> (istream& is,
            Box&     b)
{
    is >> ws;
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        is.ignore(CH_IGNORE_MAX, '(');
        is >> b.smallend ;
        is >> b.bigend ;
        IntVect v;
        is >> v;
        b.btype = IndexType(v);
        CH_assert(b.btype.ok());
        is.ignore(CH_IGNORE_MAX,')');
        b.computeBoxLen();
    }
    else if (c == '<')
    {
        is >> b.smallend;
        is >> b.bigend;
        IntVect v;
        is >> v;
        b.btype = IndexType(v);
        CH_assert(b.btype.ok());
        b.computeBoxLen();
    }
    else
        MayDay::Error("operator>>(istream&,Box&): expected \'<\'");

    if (is.fail())
        MayDay::Error("operator>>(istream&,Box&) failed");

    return is;
}

void
Box::dumpOn (ostream& strm) const
{
    strm << "Box "
         << smallend
         << " to "
         << bigend
         << " type ["
         << btype.ixType()
         << "]"
         << '\n';

    if (strm.fail())
        MayDay::Error("Box::dumpOn(ostream&) failed");
}

//
// Give smallest Box containing both Boxes.
//

Box&
Box::minBox (const Box &b)
{
  CH_assert(sameType(b));
  if (isEmpty())
  {
    *this = b;
  }
  else if (!b.isEmpty())
  {
    smallend.min(b.smallend);
    bigend.max(b.bigend);
    computeBoxLen();
  }
  return *this;
}

Box
minBox (const Box& b,
        const Box& o)
{
    CH_assert(o.sameType(b));
    if (b.isEmpty())
    {
      return(o);
    }
    else if (o.isEmpty())
    {
      return(b);
    }
    else
    {
      IntVect small = b.smallend;
      IntVect big = b.bigend;
      small.min(o.smallend);
      big.max(o.bigend);
      return Box(small,big,b.btype);
    }
}

//
// Translation functions acting on relative vectors.
//

Box
bdryBox (const Box& b,
         int        dir,
         Side::LoHiSide a_side,
         int        len)
{
  Box retval;
  if (a_side == Side::Lo)
    retval = bdryLo(b, dir, len);
  else
    retval = bdryHi(b, dir, len);

  return retval;
}

Box
bdryLo (const Box& b,
        int        dir,
        int        len)
{
  //
  // set dir'th bit to 1 = IndexType::NODE.
  //
  IndexType typ(b.btype);
  typ.set(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    int sm = low[dir];
    hi.setVal(dir,sm+len-1);
    return Box(low,hi,typ);
  }
}

Box
bdryHi (const Box& b,
        int        dir,
        int        len)
{
  //
  // Set dir'th bit to 1 = IndexType::NODE.
  //
  IndexType typ(b.btype);
  typ.set(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    unsigned int bitval = b.btype.ixType(dir);
    int bg              = hi[dir]  + 1 - bitval%2;
    low.setVal(dir,bg);
    hi.setVal(dir,bg+len-1);
    return Box(low,hi,typ);
  }
}

Box
adjCellBox (const Box& b,
            int        dir,
            Side::LoHiSide a_side,
            int        len)
{
  Box retval;
  if (a_side == Side::Lo)
    retval = adjCellLo(b, dir, len);
  else
    retval = adjCellHi(b, dir, len);

  return retval;
}

Box
adjCellLo (const Box& b,
           int        dir,
           int        len)
{
  CH_assert(len != 0);
  //
  // Set dir'th bit to 0 = IndexType::CELL.
  //
  IndexType typ(b.btype);
  typ.unset(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    int sm = low[dir];
    if (len > 0)
      {
        low.setVal(dir,sm - len);
        hi.setVal(dir,sm - 1);
      }
    else if (len < 0)
      {
        hi.setVal(dir, sm - len - 1);
      }
    return Box(low,hi,typ);
  }
}

Box
adjCellHi (const Box& b,
           int        dir,
           int        len)
{
  CH_assert(len != 0);
  //
  // Set dir'th bit to 0 = IndexType::CELL.
  //
  IndexType typ(b.btype);
  typ.unset(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    unsigned int bitval = b.btype.ixType(dir);
    int bg              = hi[dir]  + 1 - bitval%2;
    if (len > 0)
      {
        low.setVal(dir,bg);
        hi.setVal(dir,bg + len - 1);
      }
    else if (len < 0)
      {
        hi.setVal(dir,bg-1);
        low.setVal(dir, bg + len);
      }
    return Box(low,hi,typ);
  }
}

/*static*/ void Box::setTempestOutputFormat( bool b )
{
  Box::s_tempestOutputFormat = b;
}
bool Box::s_tempestOutputFormat = false;

#include "NamespaceFooter.H"
#include "UsingNamespace.H"
#include "BaseNamespaceHeader.H"

template < >
void linearIn<Box>(Box& a_outputT, const void* const a_inBuf)
{
  int* intBuf = (int*)a_inBuf;
  IntVect lo, hi;
  //output is lo[0],lo[1],...
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      lo[idir] = intBuf[idir];
      hi[idir] = intBuf[idir+SpaceDim];
    }
  if (lo <= hi)
    a_outputT = Box(lo,hi);
  else
    a_outputT = Box();
}

template < >
void linearOut<Box>(void* const a_outBuf, const Box& a_inputT)
{
  int* intBuf = (int*)a_outBuf;
  const IntVect& lo = a_inputT.smallEnd();
  const IntVect& hi = a_inputT.bigEnd();
  //output is lo[0],lo[1],...
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      intBuf[idir] = lo[idir];
      intBuf[idir+SpaceDim] = hi[idir];
    }
}

template < >
int linearSize<Box>(const Box& a_input)
{
  //box is stored as 2*spaceDim integers
  return(2*SpaceDim*sizeof(int));
}

//Vector<Box>  specialization
template < > int linearSize(const Vector<Box>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<Box>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<Box>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

//Vector<Vector<Box> >  specialization
template < > int linearSize(const Vector<Vector<Box> >& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<Vector<Box> >& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<Vector<Box> >& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

#include "BaseNamespaceFooter.H"
