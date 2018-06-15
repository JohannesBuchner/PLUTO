#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeFArrayBox.cpp
// from dMartin/Chombo-IAMR/util/FluxBox.cpp by Dan Martin, Fri, Jan 14, 2000
// petermc, 1 Nov 2000
// petermc, 26 Mar 2002 added using std::cout and using std::endl

#include "NodeFArrayBox.H"
#include "SPMD.H"
#include "NamespaceHeader.H"
using std::cout;
using std::cerr;
using std::endl;

// first do simple access functions

// ---------------------------------------------------------
const Box&
NodeFArrayBox::box() const
{
  return m_box;
}

// ---------------------------------------------------------
FArrayBox&
NodeFArrayBox::getFab()
{
  return m_fab;
}

// ---------------------------------------------------------
const FArrayBox&
NodeFArrayBox::getFab() const
{
  return m_fab;
}

// ---------------------------------------------------------
// constructors and destructors
// ---------------------------------------------------------
// ": mfab()" means initialize m_fab with default constructor.
NodeFArrayBox::NodeFArrayBox() : m_fab()
{
}

// ---------------------------------------------------------
// initialize the fab to be aliased to the input array
NodeFArrayBox::NodeFArrayBox(const Box& a_bx, int a_nComp, Real* a_alias )
  : m_fab( surroundingNodes(a_bx) ,a_nComp ,a_alias )
{
  CH_assert(a_nComp > 0);
  CH_assert(!a_bx.isEmpty());
  m_box = a_bx;
}

// ---------------------------------------------------------
NodeFArrayBox::NodeFArrayBox(const Interval& a_comps,
                             NodeFArrayBox& a_original)
  : m_fab( a_comps, a_original.m_fab )
{
  m_box = a_original.m_box;
}

// ---------------------------------------------------------
NodeFArrayBox::~NodeFArrayBox()
{

}

// ---------------------------------------------------------
// define function
void
NodeFArrayBox::define(const Box& a_bx, int a_nComp)
{
  // resize() works just like BaseFab::define() if this instance
  // hasn't already been defined
  resize(a_bx, a_nComp);
}

// ---------------------------------------------------------
// resize function
void
NodeFArrayBox::resize(const Box& a_bx, int a_nComp, Real* a_alias)
{
  CH_assert(a_nComp > 0);
  m_box = a_bx;
  Box nodeBox(surroundingNodes(a_bx));
  m_fab.resize(nodeBox, a_nComp, a_alias);
}

// ---------------------------------------------------------
void
NodeFArrayBox::copy(const NodeFArrayBox& a_src)
{
  m_fab.copy(a_src.m_fab);
}

// ---------------------------------------------------------
void
NodeFArrayBox::copy(const Box& a_regionFrom,
                    const Interval& a_Cdest,
                    const Box& a_regionTo,
                    const NodeFArrayBox& a_src,
                    const Interval& a_Csrc)
{
  CH_assert(a_regionFrom.sameSize(a_regionTo));
  IntVect translator(a_regionTo.smallEnd() - a_regionFrom.smallEnd());
  const FArrayBox& srcFab = a_src.m_fab;
  // srcFab.box() should contain surroundingNodes(a_R)
  //CH_assert(srcFab.box().contains(surroundingNodes(a_R)));
  // Box copyBox = surroundingNodes(a_R);
  Box nodesFrom = surroundingNodes(a_regionFrom);
  CH_assert(srcFab.box().contains(nodesFrom));
  Box nodesTo = surroundingNodes(a_regionTo);
  //CH_assert(m_fab.box().contains(nodesTo));
  // added by petermc, 15 Nov 2002, based on FluxBox::copy:
  // safety check -- due to node centering,
  // nodesTo may not be contained in m_fab.box().
  Box thisBox = m_fab.box();
  nodesTo &= thisBox;
  // Also need this so that nodesTo and nodesFrom are same size.
  nodesFrom &= srcFab.box(); // previously &= thisBox
  Box nodesToTranslated(nodesTo);
  nodesToTranslated.shift(-translator);
  nodesFrom &= nodesToTranslated;
  m_fab.copy(nodesFrom, a_Cdest, nodesTo, srcFab, a_Csrc);
}

// ---------------------------------------------------------
int
NodeFArrayBox::size(const Box& a_R, const Interval& a_comps) const
{
  Box nodeBox(surroundingNodes(a_R));
  return m_fab.size(nodeBox, a_comps);
}

// ---------------------------------------------------------
void
NodeFArrayBox::linearIn(void* a_buf, const Box& a_R, const Interval& a_comps)
{
  Box nodeBox(surroundingNodes(a_R));
  // added by petermc, 5 Dec 2002, to avoid problems with copyTo().
  Box thisBox = m_fab.box();
  nodeBox &= thisBox;
  m_fab.linearIn(a_buf, nodeBox, a_comps);
}

// ---------------------------------------------------------
// this is for broadcast & gather; the Box is in the msg
void
NodeFArrayBox::linearIn(const void* const a_buf)
{
  // first get the box, then the number of components, then the FAB
  char * bufptr = static_cast<char*>(const_cast<void*>(a_buf)) ;
  CH_XD::linearIn( m_box ,bufptr )  ;  bufptr += CH_XD::linearSize( m_box );
  int ncomps ;
  CH_XD::linearIn( ncomps ,bufptr ) ;  bufptr += sizeof( ncomps );
  Box nodeBox = surroundingNodes( m_box ) ;
  m_fab.resize( nodeBox ,ncomps );
  m_fab.linearIn( bufptr, nodeBox, Interval( 0,ncomps-1 ) );
}

// ---------------------------------------------------------
void
NodeFArrayBox::linearOut(void* a_buf, const Box& a_R, const Interval& a_comps) const
{
  Box nodeBox(surroundingNodes(a_R));
  // added by petermc, 5 Dec 2002, to avoid problems with copyTo().
  Box thisBox = m_fab.box();
  nodeBox &= thisBox;
  m_fab.linearOut(a_buf, nodeBox, a_comps);
}

// ---------------------------------------------------------
// this is for broadcast & gather; the Box is in the msg
void
NodeFArrayBox::linearOut(void* const a_outBuf) const
{
  char * bufptr = static_cast<char*>(a_outBuf) ;
  CH_XD::linearOut( bufptr ,m_box )          ; bufptr += CH_XD::linearSize( m_box ) ;
  CH_XD::linearOut( bufptr ,m_fab.nComp() )  ; bufptr += sizeof(int) ;
  m_fab.linearOut( bufptr ,m_fab.box() ,m_fab.interval() );
}

// ---------------------------------------------------------
int
NodeFArrayBox::linearSize( ) const
{ // the linearization contains the cell-centered Box, #components, and FAB data
  return CH_XD::linearSize(m_box) + sizeof(int) + m_fab.size( m_fab.box() ,m_fab.interval() );
}

// ---------------------------------------------------------
void
NodeFArrayBox::setVal( Real a_x )
{
    m_fab.setVal( a_x );
}

void NodeFArrayBox::setVal(Real         a_x,
                           const Box& a_bx,
                           int        a_nstart,
                           int        a_numcomp)
{
  m_fab.setVal(a_x, a_bx, a_nstart, a_numcomp);
}

#include "NamespaceFooter.H"
