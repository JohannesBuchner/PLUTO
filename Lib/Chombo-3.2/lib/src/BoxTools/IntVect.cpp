#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SPACE.H"
#include <iostream>
using std::ostream;
using std::istream;
using std::ws;

#include "MayDay.H"
#include "Misc.H"
#include "IntVect.H"
#include "parstream.H"
#include "IndexTM.H"

#include "NamespaceHeader.H"

const size_t  IntVect::IntVectSize = SpaceDim*sizeof(int);

//
// Returns IntVect which is the componentwise integer projection
// of IntVect p1 by IntVect p2.
//

ostream&
operator<< (ostream&       os,
            const IntVect& p)
{
  os << D_TERM6( '(' << p[0] , <<
                 ',' << p[1] , <<
                 ',' << p[2] , <<
                 ',' << p[3] , <<
                 ',' << p[4] , <<
                 ',' << p[5])  << ')';
    if (os.fail())
        MayDay::Error("operator<<(ostream&,IntVect&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define CH_IGNORE_MAX 100000

istream&
operator>> (istream& is,
            IntVect& p)
{
    is >> ws;
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        D_EXPR6(is.ignore(CH_IGNORE_MAX, '(') >> p[0],
                is.ignore(CH_IGNORE_MAX, ',') >> p[1],
                is.ignore(CH_IGNORE_MAX, ',') >> p[2],
                is.ignore(CH_IGNORE_MAX, '(') >> p[3],
                is.ignore(CH_IGNORE_MAX, ',') >> p[4],
                is.ignore(CH_IGNORE_MAX, ',') >> p[5]);

        is.ignore(CH_IGNORE_MAX, ')');
    }
    else if (c == '<')
    {
        D_EXPR6(is.ignore(CH_IGNORE_MAX, '<') >> p[0],
                is.ignore(CH_IGNORE_MAX, ',') >> p[1],
                is.ignore(CH_IGNORE_MAX, ',') >> p[2],
                is.ignore(CH_IGNORE_MAX, '<') >> p[3],
                is.ignore(CH_IGNORE_MAX, ',') >> p[4],
                is.ignore(CH_IGNORE_MAX, ',') >> p[5]);
        is.ignore(CH_IGNORE_MAX, '>');
    }
    else
        MayDay::Error("operator>>(istream&,IntVect&): expected \'(\' or \'<\'");

    if (is.fail())
        MayDay::Error("operator>>(istream&,IntVect&) failed");

    return is;
}

void
IntVect::printOn (ostream& os) const
{
    os << "IntVect: " << *this << '\n';
}

void
IntVect::p() const
{
    pout() << *this << '\n';
}

IntVect::IntVect(const IndexTM<int, CH_SPACEDIM>& a_tm)
{
  D_EXPR6(vect[0] = a_tm[0], vect[1] = a_tm[1], vect[2] = a_tm[2],
          vect[3] = a_tm[3], vect[4] = a_tm[4], vect[5] = a_tm[5]);
}

void
IntVect::dumpOn (ostream& os) const
{
    os << "IntVect " << *this << '\n';
}

//
// Static object initialization.
//
int IntVect::InitStatics()
{
  IntVect* pz = const_cast<IntVect*>( &IntVect::Zero );
  *pz = IntVect(D_DECL6(0,0,0,0,0,0));

  IntVect* pu = const_cast<IntVect*>( &IntVect::Unit );
  *pu = IntVect(D_DECL6(1,1,1,1,1,1));

  // No danger of IntVect::Zero and Unit not having been allocated, as ARM section
  // 3.4 says "The initialization of nonlocal static objects in a translation unit
  // is done before the first use of any function...defined in that translation
  // unit."
  //
  // Had to go through the const_cast stuff because it's nice to be able to declare
  // IntVect::Zero and IntVect::Unit as const.

  return 0; // arbitrary
}
const IntVect IntVect::Zero;
const IntVect IntVect::Unit;
static int s_dummyForIntVectCpp( IntVect::InitStatics() );
// If IntVect::Zero and IntVect::Unit were pointers, we wouldn't need this extra
// static int.  But they're objects, and the danger is that the initializations
// right above here ("IntVect IntVect::Zero;" and "IntVect IntVect::Unit;") are hit
// after the last call to IntVect::InitStatics, and in that case the
// IntVect::IntVect() constructor could redefine Zero and Unit.  In fact, the way
// things stand now, nothing bad would happen, because the IntVect::IntVect()
// constructor doesn't assign anything to any of the data members.  But we don't
// want to count on that always being the case.
#include "NamespaceFooter.H"
