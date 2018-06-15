#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// Finds sections of the ghostcell halo around a box.  Needed for
// Chombo-based implementation of Tempest's DataLayout class.
//
// This code requires the IndexTM class -- the templatized unification
// of IntVect and RealVect.
//
// The algorithm here is very inefficient (but very simple).  It should
// work for any dimension box -- 2, 3, 4, ... -- but I've only verified
// the results as correct for dimensions 2 and 3.
//
// Everything is done in terms of node-centered boxes, except for the
// Chombo interface -- the single function getBoundaryBoxes() -- which
// converts its Box argument from cell-centered to node-centered before
// passing it to the (node-centered) codimBoxes() function, and then
// converts back to cell-centered before returning.
//
// Author: Ted

#include <iostream>
#include <vector>
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "IndexTM.H"
#include "Box.H"
#include "Vector.H"
#include "parstream.H"
#include "ParmParse.H"
#include "UsingNamespace.H"

namespace CodimensionBoundary
{

/** Homegrown templatized Box class we'll need until such time as Chombo
 *  offers one.
*/
template<int DIM> class Bocks
{
  public:
    Bocks()
    {
    }

    Bocks( IndexTM<int,DIM> const & lo, IndexTM<int,DIM> const & hi )
      :
      m_smallEnd(lo),
      m_bigEnd(hi)
    {
    }

    Bocks( Bocks const & that )
      :
      m_smallEnd(that.m_smallEnd),
      m_bigEnd(that.m_bigEnd)
    {
    }

    Bocks( Box const & chomboBox )
      :
      m_smallEnd(chomboBox.smallEnd().dataPtr()),
      m_bigEnd(chomboBox.bigEnd().dataPtr())
    {
    }

    bool operator==( Bocks const & that )
    {
        return (m_smallEnd == that.m_smallEnd) && (m_bigEnd == that.m_bigEnd);
    }

    bool isDegenerate() const
    {
        for ( int d=0;d<DIM;++d ) if ( m_smallEnd[d] == m_bigEnd[d] ) return true;
        return false;
    }

    bool lexLT( Bocks<DIM> const & that ) const
    {
        if     ( this->smallEnd().lexLT( that.smallEnd() ) ) return true;
        else if ( that.smallEnd().lexLT( this->smallEnd() ) ) return false;
        else if ( this->bigEnd().lexLT( that.bigEnd() ) )     return true;
        else if ( that.bigEnd().lexLT( this->bigEnd() ) )     return false;
        else                                                 return false;
    }

    operator Box()
    {
        return Box( IntVect(smallEnd().dataPtr()), IntVect(bigEnd().dataPtr()) );
    }

    IndexTM<int,DIM> const & smallEnd() const
    {
      return m_smallEnd;
    }

    IndexTM<int,DIM> const & bigEnd() const
    {
      return m_bigEnd;
    }

    IndexTM<int,DIM> & smallEnd()
    {
      return m_smallEnd;
    }

    IndexTM<int,DIM> & bigEnd()
    {
      return m_bigEnd;
    }

  private:
    IndexTM<int,DIM> m_smallEnd;
    IndexTM<int,DIM> m_bigEnd;
};


//
// Acts as one and only point of interaction with inputs file.
//
struct InputParams
{
  InputParams( const char* a_infilename, int argc, char** argv )
    : verbose(1)
  {
      if ( argc == 0 ) return;
      ParmParse pp(0, NULL, NULL, a_infilename);
      pp.get("verbose", verbose);
  }
  int verbose;
};


template<int DIM> std::ostream&
operator<<( std::ostream & out, Bocks<DIM> const & bocks )
{
    out << "[" << bocks.smallEnd() << " " << bocks.bigEnd() << "]";
    if ( bocks.isDegenerate() ) out << " (degenerate))";
    return out;
}


template<typename T> std::ostream&
operator<<( std::ostream & out, std::vector<T> const & v )
{
    out << "{" << '\n';
    for ( unsigned i=0;i<v.size();++i )
    {
        out << v[i] << '\n';
    }
    out << "}";
    return out;
}


/** Args lo and hi are opposite corners (along the longest diagonal) of
 *  a box, but not necessarily situated on the canonical southwest-northeast
 *  diagonal.  Return an equivalent Box whose corners are at the canonical
 *  points.
*/
template<int DIM> Bocks<DIM>
canonicalizeCorners( IndexTM<int,DIM> const & lo, IndexTM<int,DIM> const & hi )
{
    IndexTM<int,DIM> goodLo, goodHi;
    for ( int i=0;i<DIM;++i )
    {
        goodLo[i] = std::min( lo[i], hi[i] );
        goodHi[i] = std::max( lo[i], hi[i] );
    }
    return Bocks<DIM>( goodLo, goodHi );
}


template<int DIM> double intvectDistance( IndexTM<int,DIM> const & iv1,
                                          IndexTM<int,DIM> const & iv2 )
{
    double ss(0.0);
    for ( int i=0;i<DIM;++i )
    {
      int diff=iv1[i] - iv2[i]; //should run faster, compiler can pipeline
      ss += diff*diff;
    }
    return pow( ss, 0.5 );
}


/** Return the box that lies between the argument boxes.
 *  The algorithm is very inefficient; it looks at all pairs of points taking
 *  one from the first box and another from the second.
 */
template<int DIM> Bocks<DIM>
interboxBox( Bocks<DIM> const & a_box1, Bocks<DIM> const & a_box2 )
{
    int box1_coords[2][DIM], box2_coords[2][DIM];
    for ( int d=0;d<DIM;++d )
    {
        box1_coords[0][d] = a_box1.smallEnd()[d];
        box1_coords[1][d] = a_box1.bigEnd()[d];
        box2_coords[0][d] = a_box2.smallEnd()[d];
        box2_coords[1][d] = a_box2.bigEnd()[d];
    }

    double smallestDistance( 1E100 );
    Bocks<DIM> result;
    int idim = 1<<DIM;  // should be same as pow(2,DIM)
    for ( int n=0; n<idim; ++n )
    {
        IndexTM<int,DIM> lo;
        for ( int d=DIM-1; d>=0; --d )
        {
            int i = (n >> d) %2;
            lo[d] = box1_coords[i][d];
        }

        for ( int m=0; m<idim; ++m )
        {
            IndexTM<int,DIM> hi;
            for ( int d=DIM-1; d>=0; --d )
            {
                int i = (m >> d) %2;
                hi[d] = box2_coords[i][d];
            }

            bool kittyCorner( true );
            for ( int d=0;d<DIM;++d )
            {
                if ( a_box1.smallEnd()[d] != a_box1.bigEnd()[d] )
                {   // Condition handles degenerate boxes when ghost=0
                    kittyCorner &= (lo[d] != hi[d]);
                }
            }
            if ( kittyCorner )
            {
                double distance( intvectDistance( lo, hi ) );
                if ( distance < smallestDistance )
                {
                    result = Bocks<DIM>( canonicalizeCorners(lo,hi) );
                    smallestDistance = distance;
                }
            }
        }
    }
    return result;
}


/** This is really codimBoxes() for the case of codim==DIM. */
template<int DIM>
std::vector< Bocks<DIM> > nodeBoxes( Bocks<DIM> const & a_box,
                                     IndexTM<int,DIM> const & a_ghosts )
{
    std::vector< Bocks<DIM> > result;

    IndexTM<int,DIM> cornerCoords[2]; // [0]==lo, [1]==hi
    cornerCoords[0] = a_box.smallEnd();
    cornerCoords[1] = a_box.bigEnd();

    int idim = 1<<DIM;  // should be same as pow(2,DIM)
    for ( int n=0; n<idim; ++n )
    {
        IndexTM<int,DIM> loCorner, hiCorner;

        for ( int d=DIM-1; d>=0; --d )
        {
            int i = (n >> d) % 2;
            loCorner[d] = cornerCoords[i][d] + (2*i - 1)*a_ghosts[d];
            hiCorner[d] = cornerCoords[i][d];
        }
        if ( loCorner.lexGT( hiCorner ) )
        {
            std::swap( loCorner, hiCorner );
        }

        result.push_back( canonicalizeCorners(loCorner,hiCorner) );
    }

    std::sort( result.begin(), result.end(), std::mem_fun_ref(&Bocks<DIM>::lexLT) );
    return result;
}


template<int DIM>
std::vector< Bocks<DIM> > codimBoxes( Bocks<DIM> const & a_box,
                                      IndexTM<int,DIM> a_ghosts,
                                      int codim )
{
    if ( DIM == codim )
    {
        return nodeBoxes( a_box, a_ghosts );
    }

    std::vector< Bocks<DIM> > result;
    std::vector< Bocks<DIM> > bookends(
        codimBoxes<DIM>(a_box,a_ghosts,codim+1));
    // If we're after edges here, then the bookends are nodes.  If we're after
    // faces, the bookends are edges.

    // Go through all unique pairs of bookends.  When they are "colinear" with
    // a box between them (i.e. when they're not at opposite corners), find
    // that box and push it onto result.

    for ( unsigned i=0; i<bookends.size(); ++i )
    {
        for ( unsigned j=i+1; j<bookends.size(); ++j )
        {
            int alignments(0);
            for ( int m=0;m<DIM;++m )
            {
                alignments +=
                  bookends[i].smallEnd()[m] == bookends[j].smallEnd()[m];
            }
            if ( ( alignments == (DIM-1) )
            &&  (!(bookends[i].isDegenerate() ^ bookends[j].isDegenerate())) )
            {
                Bocks<DIM> box( interboxBox( bookends[i], bookends[j] ) );
                if ( std::find( result.begin(), result.end(), box )
                    == result.end() )
                {
                  result.push_back( box );
                }
            }
        }
    }

    std::sort( result.begin(), result.end(), std::mem_fun_ref(&Bocks<DIM>::lexLT) );
    return result;
}

} // namespace CodimensionBoundary


//
// Chombo interface (not in namespace CodimensionBoundary anymore).
//
using namespace CodimensionBoundary;

Vector<Box> getBoundaryBoxes( const Box& a_box,
                              const IntVect& a_ghostDepth,
                              int a_codimension,
                              const InputParams& a_inputs )
{
    if ( a_inputs.verbose )
    {
      pout() << "getBoundaryBoxes( " << a_box << ", " << a_ghostDepth
                << ", " << a_codimension << ")\n";
    }

    Bocks<CH_SPACEDIM> a_bocks( a_box );
    a_bocks.bigEnd() += 1; // Converts to node-centered.
    IndexTM<int,CH_SPACEDIM> boxLo( a_bocks.smallEnd() );
    IndexTM<int,CH_SPACEDIM> boxHi( a_bocks.bigEnd() );
    Bocks<CH_SPACEDIM> bocks( boxLo, boxHi );
    IndexTM<int,CH_SPACEDIM> ghost( a_ghostDepth.dataPtr() );
    std::vector< Bocks<CH_SPACEDIM> > resultT(
         CodimensionBoundary::codimBoxes<CH_SPACEDIM>( bocks,ghost,a_codimension ));

    Vector<Box> result;
    for ( unsigned i=0;i<resultT.size();++i )
    {
        Bocks<CH_SPACEDIM> cell_centered_bocks( resultT[i] );
        cell_centered_bocks.bigEnd() -= 1; // Converts back to cell-centered.
        Box chombo_box( cell_centered_bocks );
        result.push_back( chombo_box );
    }
    return result;
}


template<int DIM> void templatizedMain( const InputParams& a_inputs )
{
    int buf1[] =
    {
      1,1,1,1,1,1
    };
    IndexTM<int,DIM> iv1( buf1 );

    int buf2[] =
    {
      3,4,5,6,7,8
    };
    IndexTM<int,DIM> iv2( buf2 );

    Bocks<DIM> box( iv1, iv2 );
    if ( a_inputs.verbose )
    {
      pout() << "Bocks: " << box << " (node-centered convention)\n";
    }

    std::vector< IndexTM<int,DIM> > ghosts_visited;

    for ( int degeneracy_order=1; degeneracy_order<=DIM; ++degeneracy_order )
    {
        for ( int g=0; g<2; ++g )  // Number of ghost cells in the x-direction.
        {
            IndexTM<int,DIM> ghosts( IndexTM<int,DIM>::Unit );
            for ( int d=0;d<degeneracy_order;++d )
            {
                ghosts[d] = g;
            }

            if ( std::find( ghosts_visited.begin(), ghosts_visited.end(), ghosts )
                != ghosts_visited.end() )
            {
                continue;
            } else
            {
                ghosts_visited.push_back( ghosts );
            }

            if ( a_inputs.verbose )
            {
              pout() << "-----------------------------------------------------\n";
              pout() << "ghosts=" << ghosts << '\n';
            }

            for ( int codim=DIM; codim>0; --codim )
            {
                std::vector< Bocks<DIM> > boundary( codimBoxes<DIM>( box, ghosts,
                                                                     codim ) );
                boundary.erase(
                    std::remove_if ( boundary.begin(), boundary.end(),
                                    std::mem_fun_ref(&Bocks<DIM>::isDegenerate) ),
                    boundary.end() );
                std::string boundary_name( "Boundary" );
                if ( codim == DIM )
                {
                    boundary_name = "Nodes";
                } else
                if ( codim == DIM-1 )
                {
                    boundary_name = "Edges";
                } else
                if ( codim == DIM-2 )
                {
                    boundary_name = "Faces";
                }
                if ( a_inputs.verbose )
                {
                  pout() << boundary_name << " (codim=" << codim << "; "
                            << boundary.size() << "):\n" << boundary << '\n';
                }
            }
        }
    }
    if ( a_inputs.verbose )
    {
      pout() << "-----------------------------------------------------\n";
      pout() << "-----------------------------------------------------\n";
    }
}


/** Chombo interface is in CodimensionBoundary.H and uses ordinary Chombo
 *  classes like Box and IntVect, instead of our templatized variants above.
*/
void testChomboInterface( const InputParams& a_inputs )
{
    IntVect iv1( D_DECL6(1,1,1,1,1,1) );
    IntVect iv2( D_DECL6(3,4,5,3,4,5) );
    Box box( iv1, iv2 );
    if ( a_inputs.verbose )
    {
      pout() << "Chombo box (cell-centered convention): " << box << '\n';
    }

    IntVect ghosts( D_DECL6(1,1,1,1,1,1) );

    for ( int codim=CH_SPACEDIM; codim>0; --codim )
    {
        Vector<Box> boundary( getBoundaryBoxes( box, ghosts, codim, a_inputs ));
        if ( a_inputs.verbose )
        {
          pout() << "Codimension " << codim << ":\n";
          for ( unsigned i=0;i<boundary.size();++i )
          {
              pout() << boundary[i] << '\n';
          }
        }
    }
}


int main( int argc, char** argv )
{
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif
    InputParams inputs( "codimensionBoundaryTest.inputs", argc, argv );

    if ( inputs.verbose ) pout() << "2D:\n";
    templatizedMain<2>( inputs );

    if ( inputs.verbose ) pout() << "\n3D:\n";
    templatizedMain<3>( inputs );

    if ( inputs.verbose ) pout() << "\n4D (unverified):\n";
    templatizedMain<4>( inputs );

    if ( inputs.verbose ) pout() << "\nChombo, hardcoded SPACEDIM=" << CH_SPACEDIM << ":\n";
    testChomboInterface( inputs );
/*
    if ( inputs.verbose ) pout() << "\n5D (unverified):\n";
    templatizedMain<5>( inputs );

    if ( inputs.verbose ) pout() << "\n6D (unverified):\n";
    templatizedMain<6>( inputs );
*/

#ifdef CH_MPI
    MPI_Finalize();
#endif
}
